$title  Aggregation Program for the GTAP9 Database

$if not set source $abort       Need to specify a source on command line: --source=xxx
$if not set target $abort       Need to specify a target on command line: --target=yyy
$if not set output $set output %target%

$set ds %source%
$include gtap9data
$include ..\defines\%target%.map

alias (ii,jj), (rr,ss);

set     gg(*)   All goods in aggregate model plus C - G - I /
                c       Household,
                g       Government consumption,
                i       Investment/;
alias (u,*);
gg(u)$ii(u) = ii(u);
abort$sum(ii$(sameas(ii,"c") or sameas(ii,"g") or sameas(ii,"i")),1) "Invalid identifier: C, G and I are reserved.";

parameters
        vom_(*,rr)      Aggretate output
        vfm_(ff,*,rr)   Endowments - Firms' purchases at market prices,
        vdfm_(ii,*,rr)  Intermediates - firms' domestic purchases at market prices,
        vifm_(ii,*,rr)  Intermediates - firms' imports at market prices,
        vxmd_(ii,rr,ss) Trade - bilateral exports at market prices,
        vst_(ii,rr)     Trade - exports for international transportation
        vtwr_(ii,jj,rr,ss)      Trade - Margins for international transportation at world prices,

        evt_(ii,rr,rr)  Volume of energy trade (mtoe),
        evd_(ii,*,rr)   Domestic energy use (mtoe),
        evi_(ii,*,rr)   Imported energy use (mtoe),
        eco2d_(ii,*,rr) CO2 emissions in domestic fuels - Mt CO2",
        eco2i_(ii,*,rr) CO2 emissions in foreign fuels - Mt CO2",

        rto_(*,rr)      Output (or income) subsidy rates
        rtf_(ff,*,rr)   Primary factor and commodity rates taxes
        rtfd_(ii,*,rr)  Firms domestic tax rates
        rtfi_(ii,*,rr)  Firms' import tax rates
        rtxs_(ii,rr,ss) Export subsidy rates
        rtms_(ii,rr,ss) Import taxes rates,

        esubd_(ii)      Elasticity of substitution (M versus D),
        esubva_(jj)     Elasticity of substitution between factors
        esubm_(ii)      Intra-import elasticity of substitution,
        etrae_(ff)      Elasticity of transformation,
        eta_(ii,rr)     Income elasticity of demand,
        epsilon_(ii,rr) Own-price elasticity of demand;


alias (i,j), (ii,jj);;
alias (r,s), (rr,ss);
set maps(s,ss), mapj(j,jj), mapg(*,*)/c.c, i.i, g.g/;
maps(r,rr) = mapr(r,rr);
mapj(j,jj) = mapi(j,jj);
mapg(i,ii) = mapi(i,ii);

file ktitle; put ktitle;

put_utility 'title' /"Aggregating vst.";
$batinclude aggr vst  i r   vst_
put_utility 'title' /"Aggregating vom.";
$batinclude aggr vom  g r   vom_
put_utility 'title' /"Aggregating vfm.";
$batinclude aggr vfm  f j r vfm_
put_utility 'title' /"Aggregating vdfm.";
$batinclude aggr vdfm i g r vdfm_
put_utility 'title' /"Aggregating vifm.";
$batinclude aggr vifm i g r vifm_
put_utility 'title' /"Aggregating vxmd.";
$batinclude aggr vxmd i r s vxmd_
put_utility 'title' /"Aggregating evd.";
$batinclude aggr evd  i g r evd_
put_utility 'title' /"Aggregating evi.";
$batinclude aggr evi  i g r evi_
put_utility 'title' /"Aggregating evt.";
$batinclude aggr evt  i r s evt_
put_utility 'title' /"Aggregating eco2d.";
$batinclude aggr eco2d  i g r eco2d_
put_utility 'title' /"Aggregating eco2i.";
$batinclude aggr eco2i  i g r eco2i_
put_utility 'title' /"Aggregating vtwr.";
$batinclude aggr vtwr i j r s vtwr_

*       First, convert tax rates into tax payments:

rto(g,r)    = rto(g,r)*vom(g,r);
rtf(f,j,r)  = rtf(f,j,r) * vfm(f,j,r);
rtfd(i,g,r) = rtfd(i,g,r) * vdfm(i,g,r);
rtfi(i,g,r) = rtfi(i,g,r) * vifm(i,g,r);
rtms(i,r,s) = rtms(i,r,s)*((1-rtxs(i,r,s)) * vxmd(i,r,s) + sum(j,vtwr(j,i,r,s)));
rtxs(i,r,s) = rtxs(i,r,s) * vxmd(i,r,s);


*       Aggregate:

put_utility 'title' /"Aggregating rto.";
$batinclude aggr rto g r   rto_
put_utility 'title' /"Aggregating rtf.";
$batinclude aggr rtf f j r rtf_
put_utility 'title' /"Aggregating rtfd.";
$batinclude aggr rtfd i g r rtfd_
put_utility 'title' /"Aggregating rtfi.";
$batinclude aggr rtfi i g r rtfi_
put_utility 'title' /"Aggregating rtxs.";
$batinclude aggr rtxs i r s rtxs_
put_utility 'title' /"Aggregating rtms.";
$batinclude aggr rtms i r s rtms_

parameter profit;
profit(gg,rr) = vom_(gg,rr)
                - sum(ii, vdfm_(ii,gg,rr) + rtfd_(ii,gg,rr))
                - sum(ii, vifm_(ii,gg,rr) + rtfi_(ii,gg,rr));

profit(gg,rr) = vom_(gg,rr) - rto_(gg,rr)
                - sum(ii, vdfm_(ii,gg,rr) + rtfd_(ii,gg,rr))
                - sum(ii, vifm_(ii,gg,rr) + rtfi_(ii,gg,rr))
                - sum(ff, vfm_(ff,gg,rr)  + rtf_(ff,gg,rr));
display profit;


*       Convert back to rates:

rto_(gg,rr)$vom_(gg,rr) = rto_(gg,rr)/vom_(gg,rr);
rtf_(ff,jj,rr)$vfm_(ff,jj,rr)  = rtf_(ff,jj,rr) / vfm_(ff,jj,rr);
rtfd_(ii,gg,rr)$ vdfm_(ii,gg,rr) = rtfd_(ii,gg,rr) / vdfm_(ii,gg,rr);
rtfi_(ii,gg,rr)$ vifm_(ii,gg,rr) = rtfi_(ii,gg,rr) / vifm_(ii,gg,rr);
rtxs_(ii,rr,ss)$ vxmd_(ii,rr,ss) = rtxs_(ii,rr,ss) / vxmd_(ii,rr,ss);
rtms_(ii,rr,ss)$((1-rtxs_(ii,rr,ss)) * vxmd_(ii,rr,ss) + sum(jj,vtwr_(jj,ii,rr,ss)))
         = rtms_(ii,rr,ss)/((1-rtxs_(ii,rr,ss)) * vxmd_(ii,rr,ss) + sum(jj,vtwr_(jj,ii,rr,ss)));

esubd_(ii)$sum(mapi(i,ii), sum((j,r), vdfm(i,j,r)+vifm(i,j,r)))
        = sum(mapi(i,ii), sum((j,r), vdfm(i,j,r)+vifm(i,j,r))*esubd(i)) /
             sum(mapi(i,ii), sum((j,r), vdfm(i,j,r)+vifm(i,j,r)));
esubva_(jj)$sum(mapi(j,jj), sum((f,r), vfm(f,j,r)))
        = sum(mapi(j,jj), sum((f,r), vfm(f,j,r)*esubva(j))) /
              sum(mapi(j,jj), sum((f,r), vfm(f,j,r)));
esubm_(ii)$sum((r,mapi(i,ii)), vim(i,r))
        = sum((r,mapi(i,ii)), vim(i,r)*esubm(i)) / sum((r,mapi(i,ii)), vim(i,r));

etrae_(ff)$sum(mapf(f,ff)$sf(f), sum((j,r), vfm(f,j,r)))
           = sum(mapf(f,ff)$sf(f), sum((j,r), vfm(f,j,r)*etrae(f))) /
             sum(mapf(f,ff)$sf(f), sum((j,r), vfm(f,j,r)));

parameter       vp(i,r) Value of private expenditure;
vp(i,r) = vdfm(i,"c",r)*(1+rtfd0(i,"c",r)) + vifm(i,"c",r)*(1+rtfi0(i,"c",r));

eta_(ii,rr)$sum((mapr(r,rr),mapi(i,ii)),vp(i,r))
        = sum((mapr(r,rr),mapi(i,ii)),vp(i,r)*eta(i,r)) /
          sum((mapr(r,rr),mapi(i,ii)),vp(i,r));
epsilon_(ii,rr)$sum((mapr(r,rr),mapi(i,ii)),vp(i,r))
        = sum((mapr(r,rr),mapi(i,ii)),vp(i,r)*epsilon(i,r)) /
          sum((mapr(r,rr),mapi(i,ii)),vp(i,r));

loop(mapf(mf,ff), etrae_(ff) = +inf;);

*---------------only for eppa-----------------------------------------------------------------------------------------------------------
parameters
SELAS_(*,layer,rr)	INITIAL SUBSTITUTION ELASTICITY MATRIX
esup_(*,rr)		fixed factor substitution elasticity
*delas			Final demand elasticity between energy and non-energy composites
*sigtrn			Elasticity between transport consumption and other consumption in final demand
d_elas_(rr)		top final demand substitution elasticity
d_elase_(rr)		substitution among energy 
d_elaso_(rr)		substitution among otherind and enerint
d_elasa_(rr)		substitution between agriculure and the rest
ESUBE_(*,rr)		Elasticity of substitution between energy inputs
sigc_			top level elasticity between energy consumption and gases
siggv_			Elasticity of substitution for ghg in vintaged sectors 
sigg_			Elasticity of substitution for ghg 
sigg0_			Elasticity of substitution for ghg 
sigu_			top level transformation elasticity between production and urban gases
pnesta_			production sectors nest a substitution elasticity   
enesta_			energy input to electricity sector nest a substitution elasticity
tnesta_			hh transport substitution elasticity between roil and the rest of own-supplied transport
tnests_			hh transport top nest elasticity (between purchased and own-supplied)
s_bc_(*,rr)		sectoral elasticity of substitution for black and organic carbon
bsigma_(rr)		Elasticity of substitution for biomass generation
boilsig_(rr)		Fixed factor elasticity for 2nd gen bio-oil
boilffg_(rr)		Fixed factor elasticity for 1st gen bio-oil
etag_(rr,*)		income elasticity of g(redefined??)
***vdfmplusvifm_(irh,ii,rr)  base year vdfm plus vifm
**price_(irh,rr,t)	Armington price with Reimer-Hertel sectors
***finalc_			HH aggregated final consumption expenditure in billion US$
**finalcp_		per capita HH aggregate final consumption expenditure in thousand US$ divided by 100
**u_(rr,t)		AIDADS utility level
***what_(irh,rr,t)	w hat
***xhat_(irh,rr,t)	x hat
***phia_(irh,rr,t)	phi in the AIDADS system
***mbs_(irh,rr,t)	marginal budget share
***elas_(irh,rr,t)	income elasticity with Reimer-Hertel sector
***elastot_(irh,rr,t)	elas converted from individual to aggregated level
***tha_(rr,t)		correction for income elasticities of other sectors
**alpha(irh)		
**beta(irh)		
**gamma(irh)		
**price07_e4_(irh,rr)  base year real price levels from EPPA4 with base year 1997
usda_(rr,*)		USDA ICP income elasticities for 2005
wsigma_(rr)		Elasticity of substitution for wind generation
bsf_(*,rr)		backstop fixed factor substitution elasticity
			
efe_(eec,*,rr)		Elasticity of substitution for ghg and energy-fe.tl (prod d)
ene_(*,rr)		Elasticity of substitution for non-energy product-ne (prod d&z&dv)
efa_(fa,eec,rr)	Elasticity of substitution for fixed factor and armington goods (prod d)
****sigma_(rr)		Elasticity of substitution for nuclear energy
neta_(rr)		parameter used to calculate sigma-nuclear elasticity
****hsigma_(rr)	Elasticity of substitution for hydroelectric
heta_(rr)		parameter used to calculate hsigma-hydroelectric elasticity	
esubi_(rr)		Elasticity of substitution for invesment 
esb_(rr)		Elasticity of substitution for household transport 
*eyt			Elasticity of substitution for international transport 
eeid_(eec,*,rr)	Elasticity of substitution for conventional intermediate energy demand 
eedf_(eec,rr)		Elasticity of substitution for conventional final energy demand  
etedf_(eec,rr)		Elasticity of substitution for conventional final energy demand -- hh transport	  
eeid_ghg_(eec,*,rr)	Elasticity of substitution for energy demand inclusive of ghg -- intermediate  
efd_ghg_(eec,*,rr)	Elasticity of substitution for energy demand inclusive of ghg -- final 
etefd_ghg_(eec,*,rr)	Elasticity of substitution for energy demand inclusive of ghg -- hh transport  

err_(rr)		Elasticity of substitution for international tansport and import-rr.tl(prod m)  
ehomm_(x,rr)		Elasticity of substitution for net imports of homogenous goods  
ehomx_(x,rr)		Elasticity of substitution for net exports of homogenous goods  
ebiom_f_(rr)		Elasticity of substitution for net imports of bio-oil  
ebiom_i_(*,rr)		Elasticity of substitution for net imports of bio-oil  
ebiox_(rr)		Elasticity of substitution for net exports of bio-oil  
eed_(rr)		Elasticity of substitution for energy and dwelling goods(prod Z) 
ew_(rr)		Elasticity of substitution for welfare  
edv_(rr)		Elasticity of substitution for dv's other input(a)   
			
ecva_(bt,rr)		Elasticity of substitution for cva (prod eb)
efva_(bt,rr)		Elasticity of substitution for fva (prod eb)
egva_(bt,rr)		Elasticity of substitution for gva (prod eb)
etdva_(bt,rr)		Elasticity of substitution for tdva (prod eb)
egtd_(bt,rr)		Elasticity of substitution for cva (prod eb)
esva_(bt,rr)		Elasticity of substitution for sva (prod eb)
eeban_(rr)		Elasticity of substitution for eb-adv-nucl 
eva_(bt,rr)		Elasticity of substitution for va (prod eb)
 			
ewiput_(rr)		Elasticity of substitution for eb-wind input-b 
ebioput_(rr)		Elasticity of substitution for eb-bioe input-b 
esiput_(rr)		Elasticity of substitution for eb-solar input-a 
ebva_(rr)		Elasticity of substitution for eb-windbio input bio cap&lab-b
egb_(rr)		Elasticity of substitution for eb-windgas input gas cap&lab-b
egc_(rr)		Elasticity of substitution for eb-windgas input for carbon-c
egf_(rr)		Elasticity of substitution for eb-windgas input for carbon&va of gas-f
eob_(rr)		Elasticity of substitution for eb-biooil input armington goods-b
**** eb still have lots parameters not included

eebv_(rr)		Elasticity of substitution for vintaged backstop production 
eeqco2_(rr)		Elasticity of substitution for transforming carbon rights to ghg rights -- regional or regional to international 
eeqghg_(ghg,rr)		Elasticity of substitution for transforming ghg rights to carbon rights -- regional or regional to international   

CH4_(rr,*)		global CH4	emissions in 2010 from edgar (Ktons)
CH4C_(rr,*)		global CH4	emissions in 2010 with combustion activity (K tons)
CH4N_(rr,*)		global CH4	emissions in 2010 with non-combustion activity (K tons)
N2O_(rr,*)		global N2O	emissions in 2010 from edgar (Ktons)
N2OC_(rr,*)		global NO2	emissions in 2010 with combustion activity (K tons)
N2ON_(rr,*)		global NO2	emissions in 2010 with non-combustion activity (K tons)
SF6_(rr,*)		global SF6      emissions in 2010 from edgar (K tons)
SF6C_(rr,*)		global SF6	emissions in 2010 with combustion activity (K tons)
SF6N_(rr,*)		global SF6	emissions in 2010 with non-combustion activity (K tons)
HFCs_(rr,*)		global HFCs     emissions in 2010 from edgar (K tons)
HFCsC_(rr,*)		global HFCs	emissions in 2010 with combustion activity (K tons)
HFCsN_(rr,*)		global HFCs	emissions in 2010 with non-combustion activity (K tons)
PFCs_(rr,*)		global PFCs     emissions in 2010 from edgar (K tons)
PFCsC_(rr,*)		global PFCs	emissions in 2010 with combustion activity (K tons)
PFCsN_(rr,*)		global PFCs	emissions in 2010 with non-combustion activity (K tons)
;

* Elasticities between non-elec and elec in cosumption, investment and 
* government.

*SCALARS  
*	NOEEC          /1.0/
*       NOEEI          /1.0/
*       NOEEG          /1.0/
*	ESUBG Elasticity in government demand /0.5/;

$ontext
**** elastisity conversion
esubva_(jj)$sum(mapi(j,jj), sum((sf,r), vfm(sf,j,r)))
        = sum(mapi(j,jj), sum((sf,r), vfm(sf,j,r)*esubva(j))) /
              sum(mapi(j,jj), sum((sf,r), vfm(sf,j,r)));
$offtext

*SELAS_(ggg,layer,rrr)	= sum((mapr(rrr,r),mapg(ggg,g)),SELAS(g,layer,r));
*SELAS_("hh",layer,rrr)	= sum(mapr(rrr,r),SELAS("hh",layer,r));
SELAS_(ii,"smm",rr)$(sum((mapr(r,rr),mapg(i,ii)),vim(i,r)))	
= sum((mapr(r,rr),mapg(i,ii)),SELAS(i,"smm",r)*vim(i,r)) /sum((mapr(r,rr),mapg(i,ii)),vim(i,r));

SELAS_(ii,"sdm",rr)$(sum((mapr(r,rr),mapi(i,ii)),sum (g,voam(i,g,r))))	
= sum((mapr(r,rr),mapi(i,ii)),SELAS(i,"sdm",r)*sum (g,voam(i,g,r))) /sum((mapr(r,rr),mapi(i,ii)),sum(g,voam(i,g,r)));

SELAS_(gg,"l_k",rr)$(sum((mapr(r,rr),mapg(g,gg)),sum (f$mf(f),vfm(f,g,r))))	
= sum((mapr(r,rr),mapg(g,gg)),SELAS(g,"l_k",r)*sum (f$mf(f),vfm(f,g,r))) /sum((mapr(r,rr),mapg(g,gg)),sum(f$mf(f),vfm(f,g,r)));

Selas_("hh","noe_el",rr)$sum((mapr(r,rr)),sum(i$(sum (mapi(i,ii)$ec(ii),1)),voam(i,"c",r)))	
= sum((mapr(r,rr)),selas("hh","noe_el",r)*sum(i$(sum (mapi(i,ii)$ec(ii),1)),voam(i,"c",r))) /sum((mapr(r,rr)),sum(i$(sum (mapi(i,ii)$ec(ii),1)),voam(i,"c",r)));

Selas_(gg,"noe_el",rr)$sum((mapr(r,rr),mapg(g,gg)),sum(i$(sum (mapi(i,ii)$ec(ii),1)),voam(i,g,r)))  
= sum((mapr(r,rr),mapg(g,gg)),selas(g,"noe_el",r)*sum(i$(sum (mapi(i,ii)$ec(ii),1)),voam(i,g,r))) /sum((mapr(r,rr),mapg(g,gg)),sum(i$(sum (mapi(i,ii)$ec(ii),1)),voam(i,g,r)));

Selas_(gg,"e_kl",rr)$(agri(gg) and sum((mapr(r,rr),mapg(g,gg)),sum(i$(sum (mapi(i,ii)$(ec(ii) or ne(ii)),1)),voam(i,g,r)))) 
= sum((mapr(r,rr),mapg(g,gg)),selas(g,"e_kl",r)*sum(i$(sum (mapi(i,ii)$(ec(ii) or ne(ii)),1)),voam(i,g,r))) /sum((mapr(r,rr),mapg(g,gg)),sum(i$(sum (mapi(i,ii)$(ec(ii) or ne(ii)),1)),voam(i,g,r)));

Selas_(gg,"e_kl",rr)$((eint(gg) or naenoe(gg)) and sum((mapr(r,rr),mapg(g,gg)),(sum(i$(sum (mapi(i,ii)$ec(ii),1)),voam(i,g,r))+sum(f$mf(f),vfm(f,g,r))))) 
= sum((mapr(r,rr),mapg(g,gg)),selas(g,"e_kl",r)*(sum(i$(sum (mapi(i,ii)$ec(ii),1)),voam(i,g,r))+sum(f$mf(f),vfm(f,g,r)))) 
/sum((mapr(r,rr),mapg(g,gg)),(sum(i$(sum (mapi(i,ii)$ec(ii),1)),voam(i,g,r))+sum(f$mf(f),vfm(f,g,r))));

*esup_(ggg,rrr)		= sum((mapr(rrr,r),mapg(ggg,g)),esup(g,r));
esup_(gg,rr)$(agri(gg) and sum((mapr(r,rr),mapg(g,gg)),(vom(g,r)-sum(mf,vfm(mf,g,r)))))	
= sum((mapr(r,rr),mapg(g,gg)),esup(g,r)*(vom(g,r)-sum(mf,vfm(mf,g,r)))) /sum((mapr(r,rr),mapg(g,gg)),(vom(g,r)-sum(mf,vfm(mf,g,r))));

esup_(gg,rr)$((eint(gg)or naenoe(gg) or enoe(gg)) and sum((mapr(r,rr),mapg(g,gg)),(vom(g,r))))	
= sum((mapr(r,rr),mapg(g,gg)),esup(g,r)*(vom(g,r))) /sum((mapr(r,rr),mapg(g,gg)),(vom(g,r)));

*d_elas_(rrr)		= sum(mapr(rrr,r),d_elas(r));
d_elas_(rr)$(sum((mapr(r,rr)),sum(i$(sum( mapi(i,ii)$nend(ii),1)),voam(i,"c",r))))	
= sum((mapr(r,rr)),d_elas(r)*sum(i$(sum (mapi(i,ii)$nend(ii),1)),voam(i,"c",r))) /sum((mapr(r,rr)),sum(i$(sum (mapi(i,ii)$nend(ii),1)),voam(i,"c",r)));
*d_elase_(rrr)		=d_elas_(rrr);
*d_elaso_(rrr)		=d_elas_(rrr);
*d_elasa_(rrr)		=d_elas_(rrr);
*ESUBE_(ggg,rrr)		= sum((mapr(rrr,r),mapg(ggg,g)),ESUBE(g,r));
ESUBE_(gg,rr)$sum((mapr(r,rr),mapg(g,gg)),sum(i$(sum (mapi(i,ii)$enoe(ii),1)),voam(i,g,r)))	
= sum((mapr(r,rr),mapg(g,gg)),ESUBE(g,r)*sum(i$(sum (mapi(i,ii)$enoe(ii),1)),voam(i,g,r))) /sum((mapr(r,rr),mapg(g,gg)),sum(i$(sum (mapi(i,ii)$enoe(ii),1)),voam(i,g,r)));

*sigc_(ghg,eec,rrr)	= sum((mapr(rrr,r),mape(eec,e)),sigc(ghg,e,r));
*sigg_(ghg,ggg,rrr)	= sum((mapr(rrr,r),mapg(ggg,g)),sigg(ghg,g,r));
*sigg_(ghg,"fd",rrr)	= sum((mapr(rrr,r)),sigg(ghg,"fd",r));
sigg_(ghg,"fd",rr)$(sum((mapr(r,rr)),vom("c",r)))	
= sum((mapr(r,rr)),sigg(ghg,"fd",r)*vom("c",r)) /sum((mapr(r,rr)),vom("c",r));

sigg_(ghg,ii,rr)$(agri(ii) and sum((mapr(r,rr),mapg(i,ii)),vom(i,r)))	
= sum((mapr(r,rr),mapg(i,ii)),sigg(ghg,i,r)*vom(i,r))/ sum((mapr(r,rr),mapg(i,ii)),vom(i,r));

sigg_("n2o",gg,rr)$(agri(gg) and sum((mapr(r,rr),mapg(g,gg)),(vom(g,r)-sum(mf,vfm(mf,g,r)))))	
= sum((mapr(r,rr),mapg(g,gg)),sigg("n2o",g,r)*(vom(g,r)-sum(mf,vfm(mf,g,r)))) /sum((mapr(r,rr),mapg(g,gg)),(vom(g,r)-sum(mf,vfm(mf,g,r))));

sigg_(ghg,gg,rr)$((eint(gg)or naenoe(gg)) and sum((mapr(r,rr),mapg(g,gg)),(vom(g,r))))	
= sum((mapr(r,rr),mapg(g,gg)),sigg(ghg,g,r)*(vom(g,r))) /sum((mapr(r,rr),mapg(g,gg)),(vom(g,r)));

sigg_(ghg,ii,rr)$((enoe(ii)) and sum((mapr(r,rr),mapg(i,ii)),(vom(i,r)-sum(sf,vfm(sf,i,r)))))	
= sum((mapr(r,rr),mapg(i,ii)),sigg(ghg,i,r)*(vom(i,r)-sum(sf,vfm(sf,i,r))))/ sum((mapr(r,rr),mapg(i,ii)),(vom(i,r)-sum(sf,vfm(sf,i,r))));


*sigg0_(ghg,uu,rrr)	=sigg(ghg,uu,rrr);
*siggv_(ggg,rrr)		= sum((mapr(rrr,r),mapg(ggg,g)),siggv(g,r));
*sigu_(ggg,rrr)		= sum((mapr(rrr,r),mapg(ggg,g)),sigu(g,r));
*sigu_("fd",rrr)		= sum((mapr(rrr,r)),sigu("fd",r));
**below only for non-electricity--------------------------------------------
sigu_(ii,rr)$(sum((mapr(r,rr),mapg(i,ii)),vom(i,r)))	
= sum((mapr(r,rr),mapg(i,ii)),sigu(i,r)*vom(i,r))/ sum((mapr(r,rr),mapg(i,ii)),vom(i,r));
*----------------------------------------------------------------------------
*s_bc_(ggg,rrr)		= sum((mapr(rrr,r),mapg(ggg,g)),s_bc(g,r));
*bsigma_(rrr)		= sum(mapr(rrr,r),bsigma(r));
*boilsig_(rrr)		= sum(mapr(rrr,r),boilsig(r));
*boilffg_(rrr)		= sum(mapr(rrr,r),boilffg(r));
*pnesta_(ggg,rrr)	= sum((mapr(rrr,r),mapg(ggg,g)),pnesta(g,r));
pnesta_(ii,rr)$(agri(ii) and sum((mapr(r,rr),mapg(i,ii)),vom(i,r)))	
= sum((mapr(r,rr),mapg(i,ii)),pnesta(i,r)*vom(i,r))/ sum((mapr(r,rr),mapg(i,ii)),vom(i,r));

pnesta_(ii,rr)$((eint(ii) or naenoe(ii) or enoe(ii)) and sum((mapr(r,rr),mapg(i,ii)),(vom(i,r)-sum(sf,vfm(sf,i,r)))))	
= sum((mapr(r,rr),mapg(i,ii)),pnesta(i,r)*(vom(i,r)-sum(sf,vfm(sf,i,r))))/ sum((mapr(r,rr),mapg(i,ii)),(vom(i,r)-sum(sf,vfm(sf,i,r))));

*pnesta_("nucl",rrr)	= sum(mapr(rrr,r),pnesta("nucl",r));

*pnesta_("hydr",rrr)	= sum(mapr(rrr,r),pnesta("hydr",r));
*enesta_(eec,rrr)	= sum((mapr(rrr,r),mape(eec,e)),enesta(e,r));
enesta_(eec,rr)$(sum((mapr(r,rr),mape(e,eec)),vom(e,r)))	
= sum((mapr(r,rr),mape(e,eec)),enesta(e,r)*vom(e,r))/sum((mapr(r,rr),mape(e,eec)),vom(e,r));

*tnesta_(rrr)		= sum(mapr(rrr,r),tnesta(r));
*tnests_(rrr)		= sum(mapr(rrr,r),tnests(r));
*etag_(rrr,ggg)		= sum((mapr(rrr,r),mapg(ggg,g)),etag(r,g));
**price_(irh,rrr,t) = sum(mapr(rrr,r),price(irh,r,t));
**finalcp_(rrr,t) = sum(mapr(rrr,r),finalcp(r,t));
**u_(rrr,t) = sum(mapr(rrr,r),u(r,t));
*usda_(rrr,iii)		= sum((mapr(rrr,r),mapi(iii,i)),usda(r,i));
*wsigma_(rrr)		= sum(mapr(rrr,r),wsigma(r));
*bsf_(uu,rrr)		= sum(mapr(rrr,r),bsf(uu,r));
			
*efe_(eec,ggg,rrr)	= sum((mapr(rrr,r),mapg(ggg,g),mape(eec,e)),efe(e,g,r));
*ene_(ggg,rrr)		= sum((mapr(rrr,r),mapg(ggg,g)),ene(g,r));
ene_(gg,rr)$sum((mapr(r,rr),mapg(g,gg)),sum(i$(sum (mapi(i,ii)$ne(ii),1)),voam(i,g,r)))		
= sum((mapr(r,rr),mapg(g,gg)),ene(g,r)*sum(i$(sum (mapi(i,ii)$ne(ii),1)),voam(i,g,r))) /sum((mapr(r,rr),mapg(g,gg)),sum(i$(sum (mapi(i,ii)$ne(ii),1)),voam(i,g,r)));

*efa_(fa,eec,rrr)	= sum((mapr(rrr,r),mape(eec,e)),efa(fa,e,r));
*neta_(rrr)		= sum(mapr(rrr,r),neta(r));

*heta_(rrr)		= sum(mapr(rrr,r),heta(r));		
*esubi_(rrr)		= sum(mapr(rrr,r),esubi(r));
esubi_(rr)$(sum(mapr(r,rr),vom("i",r)))		
= sum(mapr(r,rr),esubi(r)*vom("i",r))/sum(mapr(r,rr),vom("i",r));

*esb_(rrr)		= sum(mapr(rrr,r),esb(r));
*eeid_(eec,ggg,rrr)	= sum((mapr(rrr,r),mapg(ggg,g),mape(eec,e)),eeid(e,g,r));
*eedf_(eec,rrr)		= sum((mapr(rrr,r),mape(eec,e)),eedf(e,r));
*etedf_(eec,rrr)		= sum((mapr(rrr,r),mape(eec,e)),etedf(e,r));
*eeid_ghg_(eec,ggg,rrr)	= sum((mapr(rrr,r),mapg(ggg,g),mape(eec,e)),eeid_ghg(e,g,r));
*efd_ghg_(eec,ggg,rrr)	= sum((mapr(rrr,r),mapg(ggg,g),mape(eec,e)),efd_ghg(e,g,r));
*etefd_ghg_(eec,ggg,rrr)	= sum((mapr(rrr,r),mapg(ggg,g),mape(eec,e)),etefd_ghg(e,g,r));

*err_(rrr)		= sum(mapr(rrr,r),err(r));
err_(rr)$(sum(mapr(r,rr),(sum((i,s),vxmd(i,s,r))+sum((j,i,s),vtwr(j,i,s,r)))))		
= sum(mapr(r,rr),err(r)*(sum((i,s),vxmd(i,s,r))+sum((j,i,s),vtwr(j,i,s,r))))/ sum(mapr(r,rr),(sum((i,s),vxmd(i,s,r))+sum((j,i,s),vtwr(j,i,s,r)))) ;

*ehomm_(x,rrr)		= sum(mapr(rrr,r),ehomm(x,r));
*ehomx_(x,rrr)		= sum(mapr(rrr,r),ehomx(x,r));
*ebiom_f_(rrr)		= sum(mapr(rrr,r),ebiom_f(r));
*ebiom_i_(ggg,rrr)	= sum((mapr(rrr,r),mapg(ggg,g)),ebiom_i(g,r));
*ebiox_(rrr)		= sum(mapr(rrr,r),ebiox(r));
*eed_(rrr)		= sum(mapr(rrr,r),eed(r));
eed_(rr)$sum((mapr(r,rr)),sum(i$(sum (mapi(i,ii)$(ec(ii) or dwe(ii)),1)),voam(i,"c",r)))		
= sum((mapr(r,rr)),eed(r)*sum(i$(sum (mapi(i,ii)$(ec(ii) or dwe(ii)),1)),voam(i,"c",r))) /sum((mapr(r,rr)),sum(i$(sum (mapi(i,ii)$(ec(ii) or dwe(ii)),1)),voam(i,"c",r)));
*ew_(rrr)		= sum(mapr(rrr,r),ew(r));
ew_(rr)$(sum(mapr(r,rr),w0(r)))		
= sum(mapr(r,rr),ew(r)*w0(r))/sum(mapr(r,rr),w0(r));
*edv_(rrr)		= sum(mapr(rrr,r),edv(r));	
*ecva_(bt,rrr)		= sum(mapr(rrr,r),ecva(bt,r));
*efva_(bt,rrr)		= sum(mapr(rrr,r),efva(bt,r));
*egva_(bt,rrr)		= sum(mapr(rrr,r),egva(bt,r));
*etdva_(bt,rrr)		= sum(mapr(rrr,r),etdva(bt,r));
*egtd_(bt,rrr)		= sum(mapr(rrr,r),egtd(bt,r));
*esva_(bt,rrr)		= sum(mapr(rrr,r),esva(bt,r));
*eeban_(rrr)		= sum(mapr(rrr,r),eeban(r));	
*eva_(bt,rrr)		= sum(mapr(rrr,r),eva(bt,r)); 
*ewiput_(rrr)		= sum(mapr(rrr,r),ewiput(r)); 
*ebioput_(rrr)		= sum(mapr(rrr,r),ebioput(r));
*esiput_(rrr)		= sum(mapr(rrr,r),esiput(r));   
*ebva_(rrr)		= sum(mapr(rrr,r),ebva(r));	
*egb_(rrr)		= sum(mapr(rrr,r),egb(r)); 
*egc_(rrr)		= sum(mapr(rrr,r),egc(r)); 
*egf_(rrr)		= sum(mapr(rrr,r),egf(r));
*eob_(rrr)		= sum(mapr(rrr,r),eob(r)); 


*eebv_(rrr)		= sum(mapr(rrr,r),eebv(r));
*eeqco2_(rrr)		= sum(mapr(rrr,r),eeqco2(r));
*eeqghg_(ghg,rrr)	= sum(mapr(rrr,r),eeqghg(ghg,r));

*aggregation for ghg data

$gdxin "..\data11\GtapGhg2010.gdx"
$load CH4, CH4C, CH4N, N2O, N2OC, N2ON, SF6, SF6C, SF6N, HFCs, HFCsC, HFCsN, PFCs, PFCsC, PFCsN

put_utility 'title' /"Aggregating CH4.";
$batinclude aggr CH4  r g   CH4_

put_utility 'title' /"Aggregating CH4C.";
$batinclude aggr CH4C  r g   CH4C_

put_utility 'title' /"Aggregating CH4N.";
$batinclude aggr CH4N  r g   CH4N_

put_utility 'title' /"Aggregating N2O.";
$batinclude aggr N2O  r g   N2O_

put_utility 'title' /"Aggregating N2OC.";
$batinclude aggr N2OC  r g   N2OC_
put_utility 'title' /"Aggregating N2ON.";
$batinclude aggr N2ON  r g   N2ON_

put_utility 'title' /"Aggregating SF6.";
$batinclude aggr SF6  r g   SF6_

put_utility 'title' /"Aggregating SF6C.";
$batinclude aggr SF6C  r g   SF6C_
put_utility 'title' /"Aggregating SF6N.";
$batinclude aggr SF6N  r g   SF6N_

put_utility 'title' /"Aggregating HFCs.";
$batinclude aggr HFCs  r g   HFCs_
put_utility 'title' /"Aggregating HFCsC.";
$batinclude aggr HFCsC  r g   HFCsC_
put_utility 'title' /"Aggregating HFCsN.";
$batinclude aggr HFCsN  r g   HFCsN_

put_utility 'title' /"Aggregating PFCs.";
$batinclude aggr PFCs  r g   PFCs_
put_utility 'title' /"Aggregating PFCsC.";
$batinclude aggr PFCsC  r g   PFCsC_
put_utility 'title' /"Aggregating PFCsN.";
$batinclude aggr PFCsN  r g   PFCsN_





$if set energydata $goto energydata

put_utility 'title' /"Unloading dataset.";
execute_unload '%datadir%%output%.gdx',
        gg=g, rr=r, ff=f, ii=i,
        vfm_=vfm, vdfm_=vdfm, vifm_=vifm,vxmd_=vxmd, vst_=vst, vtwr_=vtwr,
        rto_=rto, rtf_=rtf, rtfd_=rtfd, rtfi_=rtfi, rtxs_=rtxs, rtms_=rtms,
        evd_=evd, evi_=evi, evt_=evt, eco2d_=eco2d, eco2i_=eco2i,
        esubd_=esubd, esubva_=esubva, esubm_=esubm, etrae_=etrae, eta_=eta, epsilon_=epsilon,
eec=e, x, ghg, bt, fa, layer, 
CH4_ =CH4,    
CH4C_ =CH4C,  
CH4N_ =CH4N,  
N2O_ =N2O,    
N2OC_ =N2OC,  
N2ON_ =N2ON,  
SF6_ =SF6,    
SF6C_ =SF6C,  
SF6N_ =SF6N,  
HFCs_ =HFCs,  
HFCsC_ =HFCsC,
HFCsN_ =HFCsN,
PFCs_ =PFCs,  
PFCsC_ =PFCsC,
PFCsN_ =PFCsN ,

SELAS_		=	SELAS,			  
esup_		=	esup,			  
delas,						  
sigtrn,						  
d_elas_		=	d_elas,			  
*d_elase_	=	d_elase,		  
**d_elaso_	=	d_elaso,		  
*d_elasa_	=	d_elasa,		  
ESUBE_		=	ESUBE,			  
*sigc_		=	sigc,			  
*siggv_		=	siggv,			  
sigg_		=	sigg,			  
*sigg0_		=	sigg0,			  
sigu_		=	sigu,			  
pnesta_		=	pnesta,			  
enesta_		=	enesta,			  
*tnesta_		=	tnesta,			  
*tnests_		=	tnests,			  
*s_bc_		=	s_bc,			  
*bsigma_		=	bsigma,			  
*boilsig_	=	boilsig,		  
*boilffg_	=	boilffg,		  
*etag_		=	etag,			  
 			 			  
*usda_		=	usda,			  
*wsigma_		=	wsigma,			  
*bsf_		=	bsf,			  
						  
*efe_		=	efe,			  
ene_		=	ene,			  
*efa_		=	efa,			  
*neta_		=	neta,			  
*heta_		=	heta,			        	
esubi_		=	esubi,			  
*esb_		=	esb,			  
eyt,						  
*eeid_		=	eeid,			  
*eedf_		=	eedf,			  
*etedf_		=	etedf,			  
*eeid_ghg_	=	eeid_ghg,		  
*efd_ghg_	=	efd_ghg,		  
*etefd_ghg_	=	etefd_ghg,		  
						  
err_		=	err,			  
*ehomm_		=	ehomm,			  
*ehomx_		=	ehomx,			  
*ebiom_f_	=	ebiom_f,		  
*ebiom_i_	=	ebiom_i,		  
*ebiox_		=	ebiox,			  
eed_		=	eed,			  
ew_		=	ew,			  
*edv_		=	edv,			  
						  
*ecva_		=	ecva,			  
*efva_		=	efva,			  
*egva_		=	egva,			  
*etdva_		=	etdva,			  
*egtd_		=	egtd,			  
*esva_		=	esva,			  
*eeban_		=	eeban,			  
*eva_		=	eva,			  
 						  
*ewiput_		=	ewiput,			  
*ebioput_	=	ebioput,		  
*esiput_		=	esiput,			  
*ebva_		=	ebva,			  
*egb_		=	egb, 			  
*egc_		=	egc, 			  
*egf_		=	egf,			  
*eob_		=	eob,			  
		  				  
						  
*eebv_		=	eebv,			  
*eeqco2_		=	eeqco2,			  
*eeqghg_		=	eeqghg,                               	
	
*NOEEC,		
*NOEEI,		
*NOEEG,		
ESUBG;

*---------------end of eppa----------------------------------------------------------------------------------------------------------------------------
put_utility 'title' /"All done with aggregation.";
$exit

$label energydata
put_utility 'title' /"Unloading dataset.";

parameter
        ieocarbon_(scn,ii,rr,t)         IEO carbon emissions (Mt of CO2)
        ieogdp_(scn,rr,t)               IEO gross domestic product (billion USD)
        ieoenergy_(scn,ii,*,rr,t)       IEO energy use by sector (mtoe)
        ieoele_(scn,ieo_tec, rr,t)      Power production by aggregate generaton technology
        ieocrude_(scn,rr,t)             IEO crude oil supply (mtoe)
        ieoelesup_(scn,rr,t)            IEO electricity generation and capacity (mtoe),
        unpop_(rr,t)                    UN population trajectories (millions)

        tmp1(r)                         Temporary array for aggregating one-dimensional data,
        tmp2(rr)                        Temporary array for aggregating one-dimensional data,

        tmp3(i,r)                       Temporary array for aggregating two-dimensional data,
        tmp4(ii,rr)                     Temporary array for aggregating two-dimensional data,

        tmp5(i,g,r)                     Temporary array for aggregating three-dimensional data,
        tmp6(ii,*,rr)                   Temporary array for aggregating three-dimensional data;

loop(t,
        tmp1(r) = unpop(r,t);
$batinclude aggr tmp1 r  tmp2
        unpop_(rr,t) = tmp2(rr);
);
tmp1(r) = 0; tmp2(rr) = 0;

*.$setglobal debug yes

$log Aggregating IEO Energy Projections...

loop((scn,t,ieo_tec),
        tmp1(r) = ieoele(scn,ieo_tec,r,t);
$batinclude aggr tmp1 r      tmp2
        ieoele_(scn,ieo_tec,rr,t) = tmp2(rr);
        tmp1(r) = 0; tmp2(rr) = 0;
);



loop((scn,t),
        tmp1(r) = ieogdp(scn,r,t);
$batinclude aggr tmp1 r      tmp2
        ieogdp_(scn,rr,t) = tmp2(rr);
        tmp1(r) = 0; tmp2(rr) = 0;

        tmp1(r) = ieocrude(scn,r,t);
$batinclude aggr tmp1 r      tmp2
        ieocrude_(scn,rr,t) = tmp2(rr);
        tmp1(r) = 0; tmp2(rr) = 0;

        tmp1(r) = ieoelesup(scn,r,t);
$batinclude aggr tmp1 r      tmp2
        ieoelesup_(scn,rr,t) = tmp2(rr);
        tmp1(r) = 0; tmp2(rr) = 0;

        tmp3(i,r) = ieocarbon(scn,i,r,t);
$batinclude aggr tmp3 i r      tmp4
        ieocarbon_(scn,ii,rr,t) = tmp4(ii,rr);
        tmp3(i,r) = 0; tmp4(ii,rr) = 0;

        tmp5(i,g,r) = ieoenergy(scn,i,g,r,t);
$batinclude aggr tmp5 i g r     tmp6
        ieoenergy_(scn,ii,gg,rr,t) = tmp6(ii,gg,rr);
        tmp5(i,g,r) = 0; tmp6(ii,gg,rr) = 0;
);

execute_unload '%datadir%%output%.gdx',
        gg=g, rr=r, ff=f, ii=i,
        vfm_=vfm, vdfm_=vdfm, vifm_=vifm,vxmd_=vxmd, vst_=vst, vtwr_=vtwr,
        rto_=rto, rtf_=rtf, rtfd_=rtfd, rtfi_=rtfi, rtxs_=rtxs, rtms_=rtms,
        evd_=evd, evi_=evi, evt_=evt, eco2d_=eco2d, eco2i_=eco2i,
        esubd_=esubd, esubva_=esubva, esubm_=esubm, etrae_=etrae, eta_=eta, epsilon_=epsilon,
        ieogdp_=ieogdp, ieoenergy_=ieoenergy, ieocrude_=ieocrude, ieoprice, ieocarbon_=ieocarbon, ieoelesup_=ieoelesup, ieoele_=ieoele,
        unpop_=unpop;


put_utility 'title' /"All done with aggregation.";

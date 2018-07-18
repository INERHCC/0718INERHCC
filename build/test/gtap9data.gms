$title  GAMS Code to Read a GTAP 9 Dataset

$if not set nd $set nd 0

scalar nd       Number of decimals /%nd%/;
abort$(nd<>round(nd)) "Number of decimals must be an integer";

$if not set ds $set ds gsd
$if not set yr $set yr 11
$if not set datadir $set datadir "..\data%yr%\"
$setglobal datadir %datadir%

set     g(*)    Goods plus C and G;
$if exist     %ds%.gdx $gdxin "%ds%.gdx"
$if not exist %ds%.gdx $gdxin "%datadir%%ds%.gdx"
$load g

set     i(g)    Goods
        f(*)    Factors;
$load f i
$if declared hs6 $load hs6
$if not defined r       set r(*)        Regions;
$if not defined r       $load r


set     rnum(r) Numeraire region,
        sf(f)   Sluggish primary factors (sector-specific),
        mf(f)   Mobile primary factors;

alias (r,s), (i,j);

display r;
parameters
        vfm(f,g,r)      Endowments - Firms' purchases at market prices,
        vdfm(i,g,r)     Intermediates - firms' domestic purchases at market prices,
        vifm(i,g,r)     Intermediates - firms' imports at market prices,
        vxmd(i,r,s)     Trade - bilateral exports at market prices,
        vst(i,r)        Trade - exports for international transportation
$if declared hs6        viws_hs6(i,hs6,r,s)     Bilateral trade at world prices (HS6 goods)
        vtwr(i,j,r,s)   Trade - Margins for international transportation at world prices;

$load vfm vdfm vifm vxmd vst vtwr
$if declared hs6 $load viws_hs6

if (nd>0,
        vfm(f,g,r) = vfm(f,g,r)$round(vfm(f,g,r),nd);
        vdfm(i,g,r) = vdfm(i,g,r)$round(vdfm(i,g,r),nd);
        vifm(i,g,r) = vifm(i,g,r)$round(vifm(i,g,r),nd);
        vxmd(i,r,s) = vxmd(i,r,s)$round(vxmd(i,r,s),nd);
        vst(i,r) = vst(i,r)$round(vst(i,r),nd);
        vtwr(i,j,r,s) = vtwr(i,j,r,s)$round(vtwr(i,j,r,s),nd);
$if declared hs6        viws_hs6(i,hs6,r,s) = viws_hs6(i,hs6,r,s)$round(viws_hs6(i,hs6,r,s),nd);
);

parameter
        evt(i,r,r)      Volume of energy trade (mtoe),
        evd(i,g,r)      Domestic energy use (mtoe),
        evi(i,g,r)      Imported energy use (mtoe),
        eco2d(i,g,r)    CO2 emissions in domestic fuels - Mt CO2",
        eco2i(i,g,r)    CO2 emissions in foreign fuels - Mt CO2";

$loaddc evd evi evt eco2d eco2i
if (nd>0,
        evt(i,r,s)  = evt(i,r,s)$round(evt(i,r,s), nd);
        evd(i,g,r) = evd(i,g,r)$round(evd(i,g,r), nd);
        evi(i,g,r) = evi(i,g,r)$round(evi(i,g,r), nd);
        eco2d(i,g,r) = eco2d(i,g,r)$round(eco2d(i,g,r), nd);
        eco2i(i,g,r) = eco2i(i,g,r)$round(eco2i(i,g,r), nd);
);

parameter
        rto(g,r)        Output (or income) subsidy rates
        rtf(f,g,r)      Primary factor and commodity rates taxes
        rtfd(i,g,r)     Firms domestic tax rates
        rtfi(i,g,r)     Firms' import tax rates
        rtxs(i,r,s)     Export subsidy rates
$if declared hs6        rtms_hs6(hs6,r,s)       Bilateral tariff rates
        rtms(i,r,s)     Import taxes rates;

$load rto rtf rtfd rtfi rtxs rtms
$if declared hs6 $load rtms_hs6

if (nd>0,
        rto(g,r) = rto(g,r)$round(rto(g,r),nd);
        rtf(f,g,r) = rtf(f,g,r)$round(rtf(f,g,r),nd);
        rtfd(i,g,r) = rtfd(i,g,r)$round(rtfd(i,g,r),nd);
        rtfi(i,g,r) = rtfi(i,g,r)$round(rtfi(i,g,r),nd);
        rtxs(i,r,s) = rtxs(i,r,s)$round(rtxs(i,r,s),nd);
        rtms(i,r,s) = rtms(i,r,s)$round(rtms(i,r,s),nd);
$if declared hs6  rtms_hs6(hs6,r,s) = rtms_hs6(hs6,r,s)$round(rtms_hs6(hs6,r,s),nd);
);
parameter
        esubd(i)        Elasticity of substitution (M versus D),
        esubva(g)       Elasticity of substitution between factors
        esubm(i)        Intra-import elasticity of substitution,
        etrae(f)        Elasticity of transformation,
        eta(i,r)        Income elasticity of demand,
        epsilon(i,r)    Own-price elasticity of demand;

$load esubd esubva esubm etrae eta epsilon

*       Declare some intermediate arrays which are required to
*       evaluate tax rates:

parameter       vdm(g,r)        Aggregate demand for domestic output,
                vom(g,r)        Total supply at market prices;

vdm(i,r) = sum(g, vdfm(i,g,r));
vom(i,r) = vdm(i,r) + sum(s, vxmd(i,r,s)) + vst(i,r);

parameter
        rtf0(f,g,r)     Primary factor and commodity rates taxes
        rtfd0(i,g,r)    Firms domestic tax rates
        rtfi0(i,g,r)    Firms' import tax rates
        rtxs0(i,r,s)    Export subsidy rates
        rtms0(i,r,s)    Import taxes rates;

rtf0(f,g,r) = rtf(f,g,r);
rtfd0(i,g,r) = rtfd(i,g,r);
rtfi0(i,g,r) = rtfi(i,g,r);
rtxs0(i,r,s) = rtxs(i,r,s);
rtms0(i,r,s) = rtms(i,r,s);

parameter       pvxmd(i,s,r)    Import price (power of benchmark tariff)
                pvtwr(i,s,r)    Import price for transport services;

pvxmd(i,s,r) = (1+rtms0(i,s,r)) * (1-rtxs0(i,s,r));
pvtwr(i,s,r) = 1+rtms0(i,s,r);

parameter
        vtw(j)          Aggregate international transportation services,
        vpm(r)          Aggregate private demand,
        vgm(r)          Aggregate public demand,
        vim(i,r)        Aggregate imports,
        evom(f,r)       Aggregate factor endowment at market prices,
        vb(*)           Current account balance;

vtw(j) = sum(r, vst(j,r));
vom("c",r) = sum(i, vdfm(i,"c",r)*(1+rtfd0(i,"c",r)) + vifm(i,"c",r)*(1+rtfi0(i,"c",r)))/(1-rto("c",r));
vom("g",r) = sum(i, vdfm(i,"g",r)*(1+rtfd0(i,"g",r)) + vifm(i,"g",r)*(1+rtfi0(i,"g",r)))/(1-rto("g",r));
vom("i",r) = sum(i, vdfm(i,"i",r)*(1+rtfd0(i,"i",r)) + vifm(i,"i",r)*(1+rtfi0(i,"i",r)))/(1-rto("i",r));

vdm("c",r) = vom("c",r);
vdm("g",r) = vom("g",r);
vim(i,r) =  sum(g, vifm(i,g,r));
evom(f,r) = sum(g, vfm(f,g,r));
vb(r) = vom("c",r) + vom("g",r) + vom("i",r)
        - sum(f, evom(f,r))
        - sum(g,  vom(g,r)*rto(g,r))
        - sum(g,  sum(i, vdfm(i,g,r)*rtfd(i,g,r) + vifm(i,g,r)*rtfi(i,g,r)))
        - sum(g,  sum(f, vfm(f,g,r)*rtf(f,g,r)))
        - sum((i,s), rtms(i,s,r) *  (vxmd(i,s,r) * (1-rtxs(i,s,r)) + sum(j,vtwr(j,i,s,r))))
        + sum((i,s), rtxs(i,r,s) * vxmd(i,r,s));

vb("chksum") = sum(r, vb(r));
display vb;

*       Determine which factors are sector-specific

mf(f) = yes$(1/etrae(f)=0);
sf(f) = yes$(1/etrae(f)>0);
display mf,sf;

parameter       mprofit Zero profit for m,
                yprofit Zero profit for y;

mprofit(i,r) = vim(i,r) - sum(s, pvxmd(i,s,r)*vxmd(i,s,r)+sum(j, vtwr(j,i,s,r))*pvtwr(i,s,r));
mprofit(i,r) = round(mprofit(i,r),5);
display mprofit;

yprofit(g,r) = vom(g,r)*(1-rto(g,r))
        - sum(i, vdfm(i,g,r)*(1+rtfd0(i,g,r))
               + vifm(i,g,r)*(1+rtfi0(i,g,r)))
        - sum(f, vfm(f,g,r)*(1+rtf0(f,g,r)));

yprofit(g,r) = round(yprofit(g,r),6)
display yprofit;

*       Define a numeraire region for denominating international
*       transfers:

rnum(r) = yes$(vom("c",r)=smax(s,vom("c",s)));
display rnum;
*---------only for eppa----------------------------------------------------------------------------
$gdxin "%datadir%elasticity.gdx"


set     ghg	GREEN HOUSE GAS POLLUTANTS
        bt	BACKSTOP TECHNOLOGIES
	fa	for aggregate fixed factor and armington goods (prod d)
	layer	CES NESTING HIERARCHY
	x	Homogenous goods
	e(i)	ENERGY COMMODITIES;

$load ghg bt fa layer x e

parameter
SELAS(*,layer,r)	INITIAL SUBSTITUTION ELASTICITY MATRIX
esup(*,r)		fixed factor substitution elasticity
delas			Final demand elasticity between energy and non-energy composites
sigtrn			Elasticity between transport consumption and other consumption in final demand
d_elas(r)		top final demand substitution elasticity
d_elase(r)		substitution among energy 
d_elaso(r)		substitution among otherind and enerint
d_elasa(r)		substitution between agriculure and the rest
ESUBE(*,r)		Elasticity of substitution between energy inputs
sigc			top level elasticity between energy consumption and gases
siggv			Elasticity of substitution for ghg in vintaged sectors 
sigg			Elasticity of substitution for ghg 
sigg0			Elasticity of substitution for ghg 
sigu			top level transformation elasticity between production and urban gases
pnesta			production sectors nest a substitution elasticity   
enesta			energy input to electricity sector nest a substitution elasticity
tnesta			hh transport substitution elasticity between roil and the rest of own-supplied transport
tnests			hh transport top nest elasticity (between purchased and own-supplied)
s_bc(*,r)		sectoral elasticity of substitution for black and organic carbon
bsigma(r)		Elasticity of substitution for biomass generation
boilsig(r)		Fixed factor elasticity for 2nd gen bio-oil
boilffg(r)		Fixed factor elasticity for 1st gen bio-oil
etag(r,*)		income elasticity of g(redefined??)
***vdfmplusvifm_(irh,iii,rrr)  base year vdfm plus vifm
**price_(irh,rrr,t)	Armington price with Reimer-Hertel sectors
***finalc_			HH aggregated final consumption expenditure in billion US$
**finalcp_		per capita HH aggregate final consumption expenditure in thousand US$ divided by 100
**u_(rrr,t)		AIDADS utility level
***what_(irh,rrr,t)	w hat
***xhat_(irh,rrr,t)	x hat
***phia_(irh,rrr,t)	phi in the AIDADS system
***mbs_(irh,rrr,t)	marginal budget share
***elas_(irh,rrr,t)	income elasticity with Reimer-Hertel sector
***elastot_(irh,rrr,t)	elas converted from individual to aggregated level
***tha_(rrr,t)		correction for income elasticities of other sectors
**alpha(irh)		
**beta(irh)		
**gamma(irh)		
**price07_e4_(irh,rrr)  base year real price levels from EPPA4 with base year 1997
usda(r,*)		USDA ICP income elasticities for 2005
wsigma(r)		Elasticity of substitution for wind generation
bsf(*,r)		backstop fixed factor substitution elasticity
			
efe(e,*,r)		Elasticity of substitution for ghg and energy-fe.tl (prod d)
ene(*,r)		Elasticity of substitution for non-energy product-ne (prod d&z&dv)
efa(fa,e,r)		Elasticity of substitution for fixed factor and armington goods (prod d)
****sigma(r)		Elasticity of substitution for nuclear energy
neta(r)			parameter used to calculate sigma-nuclear elasticity
****hsigma(r)		Elasticity of substitution for hydroelectric
heta(r)			parameter used to calculate hsigma-hydroelectric elasticity	
esubi(r)		Elasticity of substitution for invesment 
esb(r)			Elasticity of substitution for household transport 
eyt			Elasticity of substitution for international transport 
eeid(e,*,r)		Elasticity of substitution for conventional intermediate energy demand 
eedf(e,r)		Elasticity of substitution for conventional final energy demand  
etedf(e,r)		Elasticity of substitution for conventional final energy demand -- hh transport	  
eeid_ghg(e,*,r)		Elasticity of substitution for energy demand inclusive of ghg -- intermediate  
efd_ghg(e,*,r)		Elasticity of substitution for energy demand inclusive of ghg -- final 
etefd_ghg(e,*,r)	Elasticity of substitution for energy demand inclusive of ghg -- hh transport  

err(r)			Elasticity of substitution for international tansport and import-rr.tl(prod m)  
ehomm(x,r)		Elasticity of substitution for net imports of homogenous goods  
ehomx(x,r)		Elasticity of substitution for net exports of homogenous goods  
ebiom_f(r)		Elasticity of substitution for net imports of bio-oil  
ebiom_i(*,r)		Elasticity of substitution for net imports of bio-oil  
ebiox(r)		Elasticity of substitution for net exports of bio-oil  
eed(r)			Elasticity of substitution for energy and dwelling goods(prod Z) 
ew(r)			Elasticity of substitution for welfare  
edv(r)			Elasticity of substitution for dv's other input(a)   
			
ecva(bt,r)		Elasticity of substitution for cva (prod eb)
efva(bt,r)		Elasticity of substitution for fva (prod eb)
egva(bt,r)		Elasticity of substitution for gva (prod eb)
etdva(bt,r)		Elasticity of substitution for tdva (prod eb)
egtd(bt,r)		Elasticity of substitution for cva (prod eb)
esva(bt,r)		Elasticity of substitution for sva (prod eb)
eeban(r)		Elasticity of substitution for eb-adv-nucl 
eva(bt,r)		Elasticity of substitution for va (prod eb)
 			
ewiput(r)		Elasticity of substitution for eb-wind input-b 
ebioput(r)		Elasticity of substitution for eb-bioe input-b 
esiput(r)		Elasticity of substitution for eb-solar input-a 
ebva(r)			Elasticity of substitution for eb-windbio input bio cap&lab-b
egb(r)			Elasticity of substitution for eb-windgas input gas cap&lab-b
egc(r)			Elasticity of substitution for eb-windgas input for carbon-c
egf(r)			Elasticity of substitution for eb-windgas input for carbon&va of gas-f
eob(r)			Elasticity of substitution for eb-biooil input armington goods-b
**** eb still have lots parameters not included

eebv(r)			Elasticity of substitution for vintaged backstop production 
eeqco2(r)		Elasticity of substitution for transforming carbon rights to ghg rights -- regional or regional to international 
eeqghg(ghg,r)		Elasticity of substitution for transforming ghg rights to carbon rights -- regional or regional to international   
;

* Elasticities between non-elec and elec in cosumption, investment and 
* government.

SCALARS  
	NOEEC          
	NOEEI          
	NOEEG          
	ESUBG Elasticity in government demand ;

* add ghg data from edgar
parameters
CH4	global CH4	emissions in 2010 (K tons)
CH4C	global CH4	emissions in 2010 with combustion activity (K tons)
CH4N	global CH4	emissions in 2010 with non-combustion activity (K tons)
N2O	global N2O	emissions in 2010 (K tons)
N2OC	global NO2	emissions in 2010 with combustion activity (K tons)
N2ON	global NO2	emissions in 2010 with non-combustion activity (K tons)
SF6	global SF6      emissions in 2010 (K tons)
SF6C	global SF6	emissions in 2010 with combustion activity (K tons)
SF6N	global SF6	emissions in 2010 with non-combustion activity (K tons)
HFCs	global HFCs     emissions in 2010 (K tons)
HFCsC	global HFCs	emissions in 2010 with combustion activity (K tons)
HFCsN	global HFCs	emissions in 2010 with non-combustion activity (K tons)
PFCs	global PFCs     emissions in 2010 (K tons)
PFCsC	global PFCs	emissions in 2010 with combustion activity (K tons)
PFCsN	global PFCs	emissions in 2010 with non-combustion activity (K tons)
;					  


$load SELAS, esup,delas,sigtrn,	d_elas,	d_elase,d_elaso,d_elasa,ESUBE,sigc,siggv,sigg,
$load sigg0,sigu,pnesta,enesta,	tnesta,	tnests,	s_bc,bsigma,boilsig,boilffg,etag,usda,wsigma,bsf,	
$load efe,ene,efa,neta,heta,esubi,esb,eyt,eeid,eedf,etedf,eeid_ghg,
$load efd_ghg,etefd_ghg,err,ehomm,ehomx,ebiom_f,ebiom_i,ebiox,eed,ew,edv,ecva,efva,egva,etdva,egtd,	
$load esva,eeban,eva,ewiput,ebioput,esiput,ebva,egb,egc, 	
$load egf,eob,eebv, eeqco2,  eeqghg,NOEEC,NOEEI,NOEEG,ESUBG,

*$gdxin "..\data11\GtapGhg2010.gdx"
*$load CH4, CH4C, CH4N, N2O, N2OC, N2ON, SF6, SF6C, SF6N, HFCs, HFCsC, HFCsN, PFCs, PFCsC, PFCsN

set con  subset-consumption sector /c/
*set con(g)  subset-consumption sector /c/
*    ec(i)   Energy goods
    ec   Energy goods
                /
                coal    Coal
                oil     Crude oil
                gas     Natural gas
                roil    Refined oil products
                elec    Electricity
                /;
set	GI goverment and investment		   /g, i/;
SET     NE NON-ENERGY COMMODITIES                  /crop, live, fors, food, EINT, OTHR, serv, tran, dwe/;
set     nend non-energy commodities excluding dwe /crop, live, fors, food, EINT, OTHR, serv, tran/;
set     dwe  ownership of dwellings               /dwe/;


SET     ENOE NON-ELEC ENERGY COMMODITIES           /COAL, OIL, GAS, ROIL/;
SET     ELEC ELEC ENERGY COMMODITIES               /ELEC/;
*SET	KSI(I) CAPITAL-SPECIFIC INDUSTRIES /COAL, OIL, GAS, ROIL, ELEC, EINT/;
*SET     X(I) Homogenous goods /OIL/;
*SET     X(I) Homogenous goods /OIL, GAS/;
*SET	CGD(I) /CGD/;
*SET     OIL(I) OIL /OIL/;
*SET     GAS(I) GAS /GAS/;
*set	roil(i) /roil/;
*set	coal(i) /coal/;
*set	oil_col(enoe) /roil, coal/;
*set	refo(e) /roil/;
SET     AGRI Agriculture /crop, live, fors/,
	EINT Energy intensive /eint/;
*SET	TRANSFUEL(I)	SECTOR USING FUEL AS TRASNPORTATION /CROP, SERV, TRAN/;
SET     AENOE /COAL, OIL, GAS, ROIL, crop, live, fors, ELEC, EINT/;
SET     NAENOE /FOOD, OTHR, serv, tran, dwe/;
*SET	VEHIFUELP(E) /ROIL/;
*SET		OTE(NE) /EINT, OTHR/;
*SET		AGDEV(R) Agri-consumer developing countries /CHN, IND, REA, BRA, AFR/;

* VJK (LUC)

* LUC: Land use categories


parameter	esub(g)		Top-level elasticity indemand /c 1/
		voam(i,g,r)	value of armington goods at biginning
		rtfa0(i,g,r)    Base year effective tax rate of Armington good i use by g
		w0(r)		utility at biginning		   
		rto0(g,r)	Commodity rates taxes at biginning
		epslond(i,g,r)  co2-energy transform coeffiecient(Gton-Mtoe)
		epsloni(i,g,r)  co2-energy transform coeffiecient
		epslon(i,g,r)	co2-energy transform coeffiecient
		eco2(i,g,r)	the total co2 due to sector g using good i
		co2tax(i,g,r)	co2 tax
		clim(r)		CO2 emissions limit (Gton) 
		sclim(r)	co2 policy switch;
scalar		scale		unit chage(from Mton to Gton) /1e-3/;

rto0(g,r)    = rto(g,r);
*voam(i,g,r)  = vdfm(i,g,r) * (1+rtfd0(i,g,r)) + vifm(i,g,r) * (1+rtfi0(i,g,r));
voam(i,g,r)  = vdfm(i,g,r)  + vifm(i,g,r) ;
rtfa0(i,g,r)$voam(i,g,r)= ((vdfm(i,g,r)*(1+rtfd0(i,g,r))+vifm(i,g,r)*(1+rtfi0(i,g,r)))/voam(i,g,r))-1;

w0(r)        = vom("c",r)+vom("i",r);
*eco2d(i,g,r)$(evd(i,g,r)=0) = 0;
*eco2i(i,g,r)$(evi(i,g,r)=0) = 0;  
*eco2d(i,g,r) = eco2d(i,g,r)*scale;
*eco2i(i,g,r) = eco2i(i,g,r)*scale;
epslond(i,g,r)$evd(i,g,r) = eco2d(i,g,r) / evd(i,g,r);
epsloni(i,g,r)$evi(i,g,r) = eco2i(i,g,r) / evi(i,g,r);
*eco2(i,g,r)  = eco2d(i,g,r)+eco2i(i,g,r);
*epslon(i,g,r)$(evi(i,g,r)+evd(i,g,r)) = eco2(i,g,r) / (evi(i,g,r)+evd(i,g,r));
*epslon("coal",g,r)			= 0.1*44/12*0.24686/23.88;
*epslon("roil",g,r)			= 0.1*44/12*0.199/23.88  ;
*epslon("gas",g,r)			= 0.1*44/12*0.137/23.88  ;
co2tax(i,g,r)=0;
*clim(r)=sum((i,g),eco2d(i,g,r)+eco2i(i,g,r));
*clim(r)=sum((i,g),(evi(i,g,r)+evd(i,g,r))*epslon(i,g,r));
*eco2(i,g,r)  = (evi(i,g,r)+evd(i,g,r))*epslon(i,G,R);
sclim(r)=no;








*---------end  of  extra setting for gtapaggregation of eppa----------------------------------------------------------------------------------------
$if not set energydata $exit

*       Read the IEO energy statistics and baseline growth path:

*.$if not defined ieoscn set ieoscn     IEO scenarios    /ref,high_oil,low_oil,high_gdp,low_gdp/;

set     t       Projected and historical time periods / 2004*2050, 2055,2060,2065,2070, 2075, 2080, 2085, 2090, 2095, 2100 /,

        scn     IEO scenarios    /
                ref             Reference case,
                high_oil        High oil price,
                low_oil         Low oil price,
                high_gdp        High GDP growth,
                low_gdp         Low GDP growth /,

*.      eg(i)   Final energy goods      /ele,col,oil,gdt,gas/,
        eg(i)   Final energy goods      /ele,col,oil,gas/,

        fd      Final demand sectors    /Residential, Commercial, Industrial, ElectricPower, Transportation/,

        ieo_tec IEO technologies /
                        capacity        Generating Capacity
                        oilcap          Liquids-Fired Generating Capacity
                        gascap          Natural-Gas-Fired Generating Capacity
                        colcap          Coal-Fired Generating Capacity
                        nuccap          Nuclear Generating Capacity
                        hydrocap        Hydroelectric Renewable Generating Capacity,
                        windcap         Wind-Powered Generating Capacity,
                        geocap          Geothermal Generating Capacity,
                        othrnwcap       Other Renewable Generating Capacity,
                        solcap          Solar Generating Capacity,

                        generation      Net Electricity Generation,
                        oilgen          Net Liquids-Fired Electricity Generation,
                        gasgen          Net Natural-Gas-Fired Electricity Generation,
                        colgen          Net Coal-Fired Electricity Generation,
                        nucgen          Net Nuclear Electricity Generation,
                        hydrogen        Net Hydroelectric Generation,
                        windgen         Net Wind-Powered Electricity Generation,
                        geogen          Net Geothermal Electricity Generation,
                        othrnwgen       Net Other Renewable Electricity Generation,
                        solgen          Net Solar Electricity Generation /;


parameter
        ieocarbon(scn,*,r,t)            IEO carbon emissions by scenario (index -- %byr%=1),
        ieogdp(scn,r,t)                 IEO GDP by scenario (index -- %byr%=1),
        ieocrude(scn,r,t)               IEO crude oil supply (index -- %byr%=1),
        ieoelesup(scn,r,t)              IEO electricity supply (index -- %byr%=1),
        ieoenergy(scn,i,g,r,t)          IEO energy use by sector (index -- %byr%=1),
        ieoprice(scn,t)                 IEO oil price (index -- %byr%=1),
        ieoele(scn,ieo_tec,r,t)         IEO electricity generation and capacity (index -- %byr%=1),
        unpop(r,*)                      UN population trajectories (in millions);

$load ieocarbon ieogdp ieocrude ieoelesup ieoprice ieoenergy ieoele unpop


MINPAR=Table[{},{2}];


(* CMSSM input parameters *)

MINPAR[[1]]={{1,m0},
             {2,m12},
             {3,TanBeta},
             {4,SignumMu},
             {5,Azero}};
             
(* GMSB input parameters *)             
             
MINPAR[[2]]={{1,xiInput},   
             {2,MessengerScale},   
             {3,TanBeta},
             {4,SignumMu},
             {6,cGrav},
             {7,n5plets},
             {8,n10plets},
             {9,hsarah}
             };

UseParameterAsGUTscale = {MessengerScale};         
RealParameters = {TanBeta, xiInput, MessengerScale, n5plets, n10plets, cGrav, m0, hsarah};
ParametersToSolveTadpoles = {\[Mu],B[\[Mu]]};

RenormalizationScaleFirstGuess = m0^2 + 4 m0^2;
RenormalizationScale = Sqrt[(mq2[3, 3] + (vu^2*conj[Yu[3, 3]]*Yu[3, 3])/2)*(mu2[3, 3] + (vu^2*conj[Yu[3, 3]]*Yu[3, 3])/2)-((vd*\[Mu]*conj[Yu[3, 3]] - vu*conj[T[Yu][3, 3]])*(vd*conj[\[Mu]]*Yu[3, 3] - vu*T[Yu][3, 3]))/2];

ConditionGUTscale = g1 == g2;

BoundaryHighScale=Table[{},{2}];
BoundarySUSYScale=Table[{},{2}];
BoundaryEWSBScale=Table[{},{2}];

SelfDefinedFunctions={
"Real(dp) Function fhA(x) \n",
"Implicit None \n",
"Real(dp),Intent(in)::x \n",
"fhA= Log((1 + x)/(1 - x))/(2.*x) \n",
"End Function fhA \n \n",
"Real(dp) Function fh1(x) \n",
"Implicit None \n",
"Real(dp),Intent(in)::x \n",
"fh1= (3*((-2 + x)*Log(1 - x) - (2 + x)*Log(1 + x)))/x**4 \n",
"End Function fh1 \n \n",
"Real(dp) Function fh2(x) \n",
"Implicit None \n",
"Real(dp),Intent(in)::x \n",
"fh2=-((-2 + x**2)*(x*atanh(x) + Log(1 - x**2)))/(6.*x**2*(-1 + x**2))& \n",
"&+(24*(-1 + x)*x**2*(1 + x) - 4*Pi**2*(-1 + x)**2*x**3*(1 + x)& \n",
"&-12*(-1 + x)*x*(1 + x)*(-10 + 7*x**2)*atanh(x) + 48*x**3*atanh(x)*Log(2.0)& \n",
"&-3*(-1 + x)*(16 + 15*x - 17*x**2 - 11*x**3 + 3*x**4 + 2*x**5)*Log(1 - x)**2& \n",
"&+48*x**5*atanh(x)*Log(2*x) - 12*x**4*Log(x)*Log(1 + x)& \n",
"&-3*(1 + x)**2*(-16 + 31*x - 10*x**2 - 5*x**3 + 4*x**4)*Log(1 + x)**2& \n",
"&+Log(1 - x)*(12*x**2*(3 + x + x**3)*Log(1 + x) - 12*x**4*Log(x*(1 + x)))& \n",
"&-6*(-1 + x)*(1 + x)*(-20 + 19*x**2)*Log(1 - x**2) + 12*x**2*(2 + x**4)*Log(x)*Log(1 - x**2)& \n",
"&+12*x**2*(1 + x)**2*(2 - 2*x + x**2)*Li2(1 - x)& \n",
"&-12*(-1 + x)**2*x**2*(2 + 2*x + x**2)*Li2(1/(1 + x))& \n",
"&-24*x**3*(1 + x**2)*Li2(-1 + 2/(1 + x)))/(72.*x**2*(-1 + x**2)**2)\n",
"End Function fh2 \n \n"
};

(* GUT conditions CMSSM *)

BoundaryHighScale[[1]]={
{T[Ye], Azero*Ye},
{T[Yd], Azero*Yd},
{T[Yu], Azero*Yu},
{mq2, DIAGONAL m0^2},
{ml2, DIAGONAL m0^2},
{md2, DIAGONAL m0^2},
{mu2, DIAGONAL m0^2},
{me2, DIAGONAL m0^2},
{mHd2, m0^2},
{mHu2, m0^2},
{MassB, m12},
{MassWB,m12},
{MassG,m12}
};

(* GUT conditions GMSB *)

BoundaryHighScale[[2]]={
{MassB, gGMSB[xiInput/MessengerScale]*(n5plets + 3*n10plets)*g1^2*xiInput/(16*Pi^2)},
{MassWB,gGMSB[xiInput/MessengerScale]*(n5plets + 3*n10plets)*g2^2*xiInput/(16*Pi^2)},
{MassG, gGMSB[xiInput/MessengerScale]*(n5plets + 3*n10plets)*g3^2*xiInput/(16*Pi^2)},
(* MSSM T-terms *)
{T[Ye],0},
{T[Ye][3,3],(-1)*3*(Ye[3,3])*(hsarah^2/(4*Pi))*xiInput/(4*Pi)*(1/(2*(xiInput/MessengerScale)))*Log[(1+(xiInput/MessengerScale))/(1-(xiInput/MessengerScale))]},
{T[Yd],0},
{T[Yu],0},
(* MSSM soft masses *)
{mHd2, fGMSB[xiInput/MessengerScale]*(n5plets + 3*n10plets)*(1.5*g2^4+0.3*g1^4)*(xiInput/(16*Pi^2))^2},
{mHu2, fGMSB[xiInput/MessengerScale]*(n5plets + 3*n10plets)*(1.5*g2^4+0.3*g1^4)*(xiInput/(16*Pi^2))^2},
{mq2, DIAGONAL fGMSB[xiInput/MessengerScale]*(n5plets+3*n10plets)*(1.5*g2^4+1/30.*g1^4+8./3.*g3^4)*(xiInput/(16*Pi^2))^2},
{ml2, DIAGONAL fGMSB[xiInput/MessengerScale]*(n5plets + 3*n10plets)*(1.5*g2^4+0.3*g1^4)*(xiInput/(16*Pi^2))^2},
{ml2[3,3], ml2[3,3]+(hsarah^2/(16*Pi^2))*(xiInput^2)*(-1/6)*(xiInput/MessengerScale)^2*fh1[xiInput/MessengerScale]+((hsarah^2/(16*Pi^2))^2)*(xiInput^2)*4*fh2[xiInput/MessengerScale]},
{md2, DIAGONAL fGMSB[xiInput/MessengerScale]*(n5plets+3*n10plets)*(2/15.*g1^4+8./3.*g3^4)*(xiInput/(16*Pi^2))^2},
{mu2, DIAGONAL fGMSB[xiInput/MessengerScale]*(n5plets+3*n10plets)*(8/15.*g1^4+8./3.*g3^4)*(xiInput/(16*Pi^2))^2},
{me2, DIAGONAL fGMSB[xiInput/MessengerScale]*(n5plets+3*n10plets)*(1.2*g1^4)*(xiInput/(16*Pi^2))^2},
{me2[3,3], me2[3,3]+(hsarah^2/(16*Pi^2))*(xiInput^2)*(-1/6)*2*(xiInput/MessengerScale)^2*fh1[xiInput/MessengerScale]+((hsarah^2/(16*Pi^2))^2)*(xiInput^2)*8*fh2[xiInput/MessengerScale]}
};


InitializationValues = {
 {\[Mu],0},
 {B[\[Mu]],0}
};


UseHiggs2LoopMSSM = True;
(* UseHiggs2LoopMSSM = False; *)

QuadruplePrecision = {};

ListDecayParticles = Automatic;
ListDecayParticles3B = Automatic;



(*----------------------------------*)
(* Information for SUSY scale input *)
(*----------------------------------*)


EXTPAR={{1,M1input},
        {2,M2input},
        {3,M3input},
        {23,Muinput},
        {24,MA2input},
        {25,TanBeta}};


ParametersToSolveTadpolesLowScaleInput = {mHd2,mHu2};


BoundaryLowScaleInput={
 {MassB, M1input},
 {MassWB, M2input},
 {MassG, M3input},
 {\[Mu], Muinput},
 {B[\[Mu]], MA2input/(TanBeta + 1/TanBeta)},
 {vd,Sqrt[4 mz2/(g1^2+g2^2)]*Cos[ArcTan[TanBeta]]},
 {vu,Sqrt[4 mz2/(g1^2+g2^2)]*Sin[ArcTan[TanBeta]]}
};

(* Example for mSugra input values *)
DefaultInputValues[1] = {m0 -> 500, m12 -> 250, TanBeta -> 10, SignumMu -> 1, Azero->0 };

(* Example for GMSB input values *)
DefaultInputValues[2] = {xiInput -> 2*10^5, MessengerScale -> 10^(10), TanBeta -> 10, SignumMu -> 1, cGrav -> 1, n5plets -> 1, n10plets ->0, m0 -> 500, hsarah-> 0.9};

NeglectLoopsInvolving = {};

IncludeFineTuning = True;
FineTuningParameters={
{m0,1/2},{m12,1/2},{Azero,1/2},{\[Mu],1/2},{B[\[Mu]],1/2}
};







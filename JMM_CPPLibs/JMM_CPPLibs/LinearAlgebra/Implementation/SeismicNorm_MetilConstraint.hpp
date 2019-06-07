//
//  SeismicNorm_MetilConstraint.hpp
//  BottleNeckSeismic3
//
//  Created by Jean-Marie Mirebeau on 15/04/2019.
//  Copyright © 2019 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef SeismicNorm_MetilConstraint_hpp
#define SeismicNorm_MetilConstraint_hpp

template<typename TC, size_t VD> auto
SeismicNorm<TC,VD>::
AD2Constraint_Metil(const VectorType & p) const -> AD2Type {
	
	const auto & p_coeffs = p;
	const auto & hooke_coeffs = hookeTensor.data;
	const int offset_in_p_coeffs = 0;
	const int offset_in_res = 0;
	
	std::array<ComponentType,10> res;
	
	ComponentType R0 = p_coeffs[ offset_in_p_coeffs + 0 ]; ComponentType R1 = hooke_coeffs[ 15 ]; ComponentType R2 = R1*R0; ComponentType R3 = hooke_coeffs[ 1 ]; ComponentType R4 = hooke_coeffs[ 20 ]; ComponentType R5 = R3+R4;
	ComponentType R6 = p_coeffs[ offset_in_p_coeffs + 1 ]; ComponentType R7 = R5*R6; ComponentType R8 = R2+R7; ComponentType R9 = R0*R8; ComponentType R10 = hooke_coeffs[ 16 ]; ComponentType R11 = R10*R6;
	ComponentType R12 = R6*R11; ComponentType R13 = p_coeffs[ offset_in_p_coeffs + 2 ]; ComponentType R14 = hooke_coeffs[ 11 ]; ComponentType R15 = hooke_coeffs[ 18 ]; ComponentType R16 = R14+R15; ComponentType R17 = R16*R6;
	ComponentType R18 = hooke_coeffs[ 6 ]; ComponentType R19 = hooke_coeffs[ 19 ]; ComponentType R20 = R18+R19; ComponentType R21 = R0*R20; ComponentType R22 = hooke_coeffs[ 13 ]; ComponentType R23 = R22*R13;
	ComponentType R24 = R17+R21+R23; ComponentType R25 = R13*R24; ComponentType R26 = R9+R12+R25; int R27 = -2; ComponentType R28 = hooke_coeffs[ 7 ]; ComponentType R29 = R28*R6;
	ComponentType R30 = R0*R16; ComponentType R31 = hooke_coeffs[ 4 ]; ComponentType R32 = hooke_coeffs[ 9 ]; ComponentType R33 = R31+R32; ComponentType R34 = R33*R13; ComponentType R35 = R30+R34;
	ComponentType R36 = R29+R35; ComponentType R37 = R6*R36; ComponentType R38 = hooke_coeffs[ 8 ]; ComponentType R39 = R38*R13; ComponentType R40 = R13*R39; ComponentType R41 = R19*R0;
	ComponentType R42 = hooke_coeffs[ 17 ]; ComponentType R43 = R22+R42; ComponentType R44 = R43*R13; ComponentType R45 = R41+R44; ComponentType R46 = R0*R45; ComponentType R47 = R37+R40+R46;
	ComponentType R48 = hooke_coeffs[ 12 ]; ComponentType R49 = R48*R13; ComponentType R50 = R13*R49; ComponentType R51 = R15*R6; ComponentType R52 = R51+R44; ComponentType R53 = R21+R52;
	ComponentType R54 = R6*R53; ComponentType R55 = hooke_coeffs[ 10 ]; ComponentType R56 = R55*R0; ComponentType R57 = hooke_coeffs[ 3 ]; ComponentType R58 = hooke_coeffs[ 14 ]; ComponentType R59 = R57+R58;
	ComponentType R60 = R59*R13; ComponentType R61 = R56+R60; ComponentType R62 = R0*R61; ComponentType R63 = R50+R54+R62; ComponentType R64 = R27*R47*R63; int R65 = -1;
	ComponentType R66 = R58*R0; ComponentType R67 = R0*R66; ComponentType R68 = hooke_coeffs[ 5 ]; ComponentType R69 = R68*R13; int R70 = 2; ComponentType R71 = R48*R0;
	ComponentType R72 = R38*R6; ComponentType R73 = R71+R72; ComponentType R74 = R70*R73; ComponentType R75 = R69+R74; ComponentType R76 = R13*R75; ComponentType R77 = R32*R6;
	ComponentType R78 = R0*R70; ComponentType R79 = R22*R78; ComponentType R80 = R77+R79; ComponentType R81 = R6*R80; ComponentType R82 = R65+R67+R76+R81; ComponentType R83 = R82*R26;
	ComponentType R84 = R64+R83; ComponentType R85 = R26*R84; ComponentType R86 = R47*R47; ComponentType R87 = R4*R0; ComponentType R88 = R0*R87; ComponentType R89 = R32*R13;
	ComponentType R90 = R15*R0; ComponentType R91 = R29+R90; ComponentType R92 = R70*R91; ComponentType R93 = R89+R92; ComponentType R94 = R13*R93; ComponentType R95 = hooke_coeffs[ 2 ];
	ComponentType R96 = R95*R6; ComponentType R97 = R10*R78; ComponentType R98 = R96+R97; ComponentType R99 = R6*R98; ComponentType R100 = R65+R88+R94+R99; ComponentType R101 = R82*R100;
	ComponentType R102 = R86-R101; ComponentType R103 = R13*R13; ComponentType R104 = R58*R103; ComponentType R105 = hooke_coeffs[ 0 ]; ComponentType R106 = R105*R0; ComponentType R107 = R1*R6;
	ComponentType R108 = R55*R13; ComponentType R109 = R107+R108; ComponentType R110 = R70*R109; ComponentType R111 = R106+R110; ComponentType R112 = R0*R111; ComponentType R113 = R4*R6;
	ComponentType R114 = R13*R70; ComponentType R115 = R19*R114; ComponentType R116 = R113+R115; ComponentType R117 = R6*R116; ComponentType R118 = R65+R104+R112+R117; ComponentType R119 = R102*R118;
	ComponentType R120 = R6*R53; ComponentType R121 = R62+R50+R120; ComponentType R122 = R100*R63; ComponentType R123 = R121*R122; ComponentType R124 = R85+R119+R123; res[ offset_in_res + 0 ] = R124;
	ComponentType R125 = R1*R78; ComponentType R126 = R20*R13; ComponentType R127 = R7+R125+R126; ComponentType R128 = R82*R127; ComponentType R129 = R19*R78; ComponentType R130 = R17+R44+R129;
	ComponentType R131 = R130*R63; ComponentType R132 = R20*R6; ComponentType R133 = R55*R78; ComponentType R134 = R132+R60+R133; ComponentType R135 = R47*R134; ComponentType R136 = R131+R135;
	ComponentType R137 = R22*R6; ComponentType R138 = R137+R49+R66; ComponentType R139 = R26*R138; ComponentType R140 = R136-R139; ComponentType R141 = R27*R140; ComponentType R142 = R128+R141;
	ComponentType R143 = R26*R142; ComponentType R144 = R127*R84; ComponentType R145 = R134*R100; ComponentType R146 = R15*R13; ComponentType R147 = R146+R11+R87; ComponentType R148 = R63*R70;
	ComponentType R149 = R147*R148; ComponentType R150 = R145+R149; ComponentType R151 = R150*R121; ComponentType R152 = R134*R122; ComponentType R153 = R82*R147; ComponentType R154 = R138*R100;
	ComponentType R155 = R153+R154; ComponentType R156 = R47*R130; ComponentType R157 = R155-R156; ComponentType R158 = R118*R157; ComponentType R159 = R106+R109; ComponentType R160 = R102*R159;
	ComponentType R161 = R158-R160; ComponentType R162 = R27*R161; ComponentType R163 = R143+R144+R151+R152+R162; res[ offset_in_res + 1 ] = R163; ComponentType R164 = R0*R5; ComponentType R165 = R70*R11;
	ComponentType R166 = R16*R13; ComponentType R167 = R164+R165+R166; ComponentType R168 = R82*R167; ComponentType R169 = R6*R70; ComponentType R170 = R28*R169; ComponentType R171 = R170+R35;
	ComponentType R172 = R171*R63; ComponentType R173 = R15*R169; ComponentType R174 = R173+R44; ComponentType R175 = R21+R174; ComponentType R176 = R47*R175; ComponentType R177 = R172+R176;
	ComponentType R178 = R22*R0; ComponentType R179 = R178+R39+R77; ComponentType R180 = R26*R179; ComponentType R181 = R177-R180; ComponentType R182 = R27*R181; ComponentType R183 = R168+R182;
	ComponentType R184 = R26*R183; ComponentType R185 = R167*R84; ComponentType R186 = R175*R100; ComponentType R187 = R10*R0; ComponentType R188 = R28*R13; ComponentType R189 = R187+R188+R96;
	ComponentType R190 = R189*R148; ComponentType R191 = R186+R190; ComponentType R192 = R191*R121; ComponentType R193 = R175*R122; ComponentType R194 = R82*R189; ComponentType R195 = R179*R100;
	ComponentType R196 = R194+R195; ComponentType R197 = R171*R47; ComponentType R198 = R196-R197; ComponentType R199 = R118*R198; ComponentType R200 = R19*R13; ComponentType R201 = R200+R2+R113;
	ComponentType R202 = R102*R201; ComponentType R203 = R199-R202; ComponentType R204 = R27*R203; ComponentType R205 = R184+R185+R192+R193+R204; res[ offset_in_res + 2 ] = R205; ComponentType R206 = R22*R114;
	ComponentType R207 = R21+R17+R206; ComponentType R208 = R82*R207; ComponentType R209 = R70*R39; ComponentType R210 = R33*R6; ComponentType R211 = R0*R43; ComponentType R212 = R209+R210+R211;
	ComponentType R213 = R212*R63; ComponentType R214 = R70*R49; ComponentType R215 = R43*R6; ComponentType R216 = R0*R59; ComponentType R217 = R214+R215+R216; ComponentType R218 = R47*R217;
	ComponentType R219 = R213+R218; ComponentType R220 = R69+R73; ComponentType R221 = R26*R220; ComponentType R222 = R219-R221; ComponentType R223 = R27*R222; ComponentType R224 = R208+R223;
	ComponentType R225 = R26*R224; ComponentType R226 = R207*R84; ComponentType R227 = R217*R100; ComponentType R228 = R90+R29+R89; ComponentType R229 = R228*R148; ComponentType R230 = R227+R229;
	ComponentType R231 = R230*R121; ComponentType R232 = R217*R122; ComponentType R233 = R82*R228; ComponentType R234 = R220*R100; ComponentType R235 = R233+R234; ComponentType R236 = R47*R212;
	ComponentType R237 = R235-R236; ComponentType R238 = R118*R237; ComponentType R239 = R19*R6; ComponentType R240 = R58*R13; ComponentType R241 = R239+R56+R240; ComponentType R242 = R102*R241;
	ComponentType R243 = R238-R242; ComponentType R244 = R27*R243; ComponentType R245 = R225+R226+R231+R232+R244; res[ offset_in_res + 3 ] = R245; int R246 = -4; ComponentType R247 = R55*R47;
	ComponentType R248 = R130*R134; ComponentType R249 = R19*R63; ComponentType R250 = R247+R248+R249; ComponentType R251 = R127*R138; ComponentType R252 = R250-R251; ComponentType R253 = R246*R252;
	ComponentType R254 = R1*R82; ComponentType R255 = R58*R26; ComponentType R256 = R254+R255; ComponentType R257 = R70*R256; ComponentType R258 = R253+R257; ComponentType R259 = R26*R258;
	int R260 = -8; ComponentType R261 = R260*R138*R147; int R262 = 4; ComponentType R263 = R47*R262; ComponentType R264 = R19*R263; ComponentType R265 = R4*R82;
	ComponentType R266 = R58*R100; ComponentType R267 = R265+R266; ComponentType R268 = R130*R130; ComponentType R269 = R267-R268; ComponentType R270 = R27*R269; ComponentType R271 = R261+R264+R270;
	ComponentType R272 = R118*R271; ComponentType R273 = R260*R157*R159; ComponentType R274 = R262*R134*R147; ComponentType R275 = R55*R100; ComponentType R276 = R4*R63; ComponentType R277 = R275+R276;
	ComponentType R278 = R70*R277; ComponentType R279 = R274+R278; ComponentType R280 = R121*R279; ComponentType R281 = R55*R122; ComponentType R282 = R127*R142; ComponentType R283 = R1*R84;
	ComponentType R284 = R105*R102; ComponentType R285 = R134*R150; ComponentType R286 = R281+R282+R283+R284+R285; ComponentType R287 = R70*R286; ComponentType R288 = R259+R272+R273+R280+R287; res[ offset_in_res + 4 ] = R288;
	ComponentType R289 = R127*R183; ComponentType R290 = R167*R142; ComponentType R291 = R5*R82; ComponentType R292 = R20*R47; ComponentType R293 = R16*R63; ComponentType R294 = R130*R175;
	ComponentType R295 = R171*R134; ComponentType R296 = R292+R293+R294+R295; ComponentType R297 = R22*R26; ComponentType R298 = R167*R138; ComponentType R299 = R127*R179; ComponentType R300 = R297+R298+R299;
	ComponentType R301 = R296-R300; ComponentType R302 = R27*R301; ComponentType R303 = R291+R302; ComponentType R304 = R26*R303; ComponentType R305 = R5*R84; ComponentType R306 = R10*R82;
	ComponentType R307 = R22*R100; ComponentType R308 = R306+R307; ComponentType R309 = R171*R130; ComponentType R310 = R16*R47; ComponentType R311 = R309+R310; ComponentType R312 = R308-R311;
	ComponentType R313 = R27*R312; ComponentType R314 = R179*R147; ComponentType R315 = R138*R189; ComponentType R316 = R314+R315; ComponentType R317 = R246*R316; ComponentType R318 = R313+R317;
	ComponentType R319 = R118*R318; ComponentType R320 = R102*R70; ComponentType R321 = R1*R320; ComponentType R322 = R150*R175; ComponentType R323 = R134*R191; ComponentType R324 = R20*R100;
	ComponentType R325 = R10*R63; ComponentType R326 = R175*R147; ComponentType R327 = R134*R189; ComponentType R328 = R325+R326+R327; ComponentType R329 = R70*R328; ComponentType R330 = R324+R329;
	ComponentType R331 = R121*R330; ComponentType R332 = R20*R122; ComponentType R333 = R157*R201; ComponentType R334 = R159*R198; ComponentType R335 = R333+R334; ComponentType R336 = R246*R335;
	ComponentType R337 = R289+R290+R304+R305+R319+R321+R322+R323+R331+R332+R336; res[ offset_in_res + 5 ] = R337; ComponentType R338 = R15*R47; ComponentType R339 = R171*R175; ComponentType R340 = R28*R63; ComponentType R341 = R338+R339+R340;
	ComponentType R342 = R167*R179; ComponentType R343 = R341-R342; ComponentType R344 = R246*R343; ComponentType R345 = R32*R26; ComponentType R346 = R306+R345; ComponentType R347 = R70*R346;
	ComponentType R348 = R344+R347; ComponentType R349 = R26*R348; ComponentType R350 = R260*R179*R189; ComponentType R351 = R28*R263; ComponentType R352 = R95*R82; ComponentType R353 = R32*R100;
	ComponentType R354 = R352+R353; ComponentType R355 = R171*R171; ComponentType R356 = R354-R355; ComponentType R357 = R27*R356; ComponentType R358 = R350+R351+R357; ComponentType R359 = R118*R358;
	ComponentType R360 = R260*R198*R201; ComponentType R361 = R262*R175*R189; ComponentType R362 = R95*R63; ComponentType R363 = R15*R100; ComponentType R364 = R362+R363; ComponentType R365 = R70*R364;
	ComponentType R366 = R361+R365; ComponentType R367 = R121*R366; ComponentType R368 = R167*R183; ComponentType R369 = R10*R84; ComponentType R370 = R15*R122; ComponentType R371 = R4*R102;
	ComponentType R372 = R191*R175; ComponentType R373 = R368+R369+R370+R371+R372; ComponentType R374 = R70*R373; ComponentType R375 = R349+R359+R360+R367+R374; res[ offset_in_res + 6 ] = R375; ComponentType R376 = R127*R224;
	ComponentType R377 = R207*R142; ComponentType R378 = R20*R82; ComponentType R379 = R59*R47; ComponentType R380 = R43*R63; ComponentType R381 = R130*R217; ComponentType R382 = R212*R134;
	ComponentType R383 = R379+R380+R381+R382; ComponentType R384 = R48*R26; ComponentType R385 = R207*R138; ComponentType R386 = R127*R220; ComponentType R387 = R384+R385+R386; ComponentType R388 = R383-R387;
	ComponentType R389 = R27*R388; ComponentType R390 = R378+R389; ComponentType R391 = R26*R390; ComponentType R392 = R20*R84; ComponentType R393 = R55*R320; ComponentType R394 = R48*R100;
	ComponentType R395 = R15*R82; ComponentType R396 = R394+R395; ComponentType R397 = R130*R212; ComponentType R398 = R43*R47; ComponentType R399 = R397+R398; ComponentType R400 = R396-R399;
	ComponentType R401 = R27*R400; ComponentType R402 = R220*R147; ComponentType R403 = R138*R228; ComponentType R404 = R402+R403; ComponentType R405 = R246*R404; ComponentType R406 = R401+R405;
	ComponentType R407 = R118*R406; ComponentType R408 = R217*R150; ComponentType R409 = R134*R230; ComponentType R410 = R59*R100; ComponentType R411 = R15*R63; ComponentType R412 = R217*R147;
	ComponentType R413 = R134*R228; ComponentType R414 = R411+R412+R413; ComponentType R415 = R70*R414; ComponentType R416 = R410+R415; ComponentType R417 = R121*R416; ComponentType R418 = R59*R122;
	ComponentType R419 = R157*R241; ComponentType R420 = R159*R237; ComponentType R421 = R419+R420; ComponentType R422 = R246*R421; ComponentType R423 = R376+R377+R391+R392+R393+R407+R408+R409+R417+R418+R422; res[ offset_in_res + 7 ] = R423;
	ComponentType R424 = R16*R84; ComponentType R425 = R16*R82; ComponentType R426 = R33*R63; ComponentType R427 = R171*R217; ComponentType R428 = R212*R175; ComponentType R429 = R398+R426+R427+R428;
	ComponentType R430 = R38*R26; ComponentType R431 = R207*R179; ComponentType R432 = R167*R220; ComponentType R433 = R430+R431+R432; ComponentType R434 = R429-R433; ComponentType R435 = R27*R434;
	ComponentType R436 = R425+R435; ComponentType R437 = R26*R436; ComponentType R438 = R167*R224; ComponentType R439 = R207*R183; ComponentType R440 = R19*R320; ComponentType R441 = R28*R82;
	ComponentType R442 = R38*R100; ComponentType R443 = R441+R442; ComponentType R444 = R171*R212; ComponentType R445 = R33*R47; ComponentType R446 = R444+R445; ComponentType R447 = R443-R446;
	ComponentType R448 = R27*R447; ComponentType R449 = R220*R189; ComponentType R450 = R179*R228; ComponentType R451 = R449+R450; ComponentType R452 = R246*R451; ComponentType R453 = R448+R452;
	ComponentType R454 = R118*R453; ComponentType R455 = R217*R191; ComponentType R456 = R230*R175; ComponentType R457 = R43*R100; ComponentType R458 = R175*R228; ComponentType R459 = R217*R189;
	ComponentType R460 = R458+R340+R459; ComponentType R461 = R70*R460; ComponentType R462 = R457+R461; ComponentType R463 = R121*R462; ComponentType R464 = R43*R122; ComponentType R465 = R198*R241;
	ComponentType R466 = R201*R237; ComponentType R467 = R465+R466; ComponentType R468 = R246*R467; ComponentType R469 = R424+R437+R438+R439+R440+R454+R455+R456+R463+R464+R468; res[ offset_in_res + 8 ] = R469; ComponentType R470 = R22*R84;
	ComponentType R471 = R48*R122; ComponentType R472 = R207*R224; ComponentType R473 = R58*R102; ComponentType R474 = R217*R230; ComponentType R475 = R470+R471+R472+R473+R474; ComponentType R476 = R70*R475;
	ComponentType R477 = R48*R47; ComponentType R478 = R212*R217; ComponentType R479 = R38*R63; ComponentType R480 = R477+R478+R479; ComponentType R481 = R207*R220; ComponentType R482 = R480-R481;
	ComponentType R483 = R246*R482; ComponentType R484 = R22*R82; ComponentType R485 = R68*R26; ComponentType R486 = R484+R485; ComponentType R487 = R70*R486; ComponentType R488 = R483+R487;
	ComponentType R489 = R26*R488; ComponentType R490 = R260*R220*R228; ComponentType R491 = R38*R263; ComponentType R492 = R32*R82; ComponentType R493 = R68*R100; ComponentType R494 = R492+R493;
	ComponentType R495 = R212*R212; ComponentType R496 = R494-R495; ComponentType R497 = R27*R496; ComponentType R498 = R490+R491+R497; ComponentType R499 = R118*R498; ComponentType R500 = R260*R237*R241;
	ComponentType R501 = R262*R217*R228; ComponentType R502 = R32*R63; ComponentType R503 = R394+R502; ComponentType R504 = R70*R503; ComponentType R505 = R501+R504; ComponentType R506 = R121*R505;
	ComponentType R507 = R476+R489+R499+R500+R506; res[ offset_in_res + 9 ] = R507;  /* 586 instructions */
	
	AD2Type result;
	auto resIt = res.begin();
	result.x = *resIt; ++resIt;
	for(int i=0; i<Dimension; ++i) {result.v[i] = *resIt; ++resIt;}
	for(int i=0; i<(Dimension*(Dimension+1))/2; ++i){result.m.data[i] = *resIt; ++resIt;}
	
	return result;
}


#endif /* SeismicNorm_MetilConstraint_hpp */

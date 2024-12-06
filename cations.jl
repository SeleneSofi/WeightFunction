Elim = -599.0175794;
ζ::Int16 = 19
ne::Int16 = 18
deltz = +2*log(19/18)
qt = [
    ETabn(-5+deltz,20,50),
    ETabn(-5+deltz,15,40)
]

deltz = +2*log(19/18)
x0 = [
    nSig(1, 1, 1.0367695364495866, 8.588072619876678e-12, Float64[], [0.001276837951486954, 0.34201177415856276, 0.6226922311756361, 0.03401915671431427], [2.7996198519362414, 4.7160286475148325, 5.182955163343529, 6.43005828877047].+deltz)
    nSig(2, 1, 1.8164287483582775, 0.22280876454319756, [2.8344173833300093].+deltz, [0.173196526504902, 0.3574556409679281, 0.4590691022249334, 0.01027873030223656], [2.362077836615896, 3.006583583628401, 3.903486390972554, 5.967283236044136].+deltz)
    nSig(3, 1, 1.9277084437825913, 0.19501595158806834, [0.41324396270039504, 2.820390624753934].+deltz, [0.0030266717530555588, 0.0533212494296452, 0.5176593839406324, 0.4259926948766669], [-0.123195973688624, 0.7034185946618937, 1.6733027359323076, 3.4262960763185886].+deltz)
    nSig(2, 2, 1.175173203207956, 0.13800245967980745, Float64[], [8.934255954285922e-5, 0.48884977099607013, 0.4513331974614988, 0.05972768898288823], [0.16953043820732777, 2.9351988058236023, 4.063902885873775, 5.681658738455484].+deltz)
    nSig(3, 2, 1.787623633012238, 0.21298255030688235, [1.0785040392562348].+deltz, [0.008185227511168926, 0.10500337321483771, 0.6215390746633012, 0.2652723246106921], [-0.45769323473291124, 0.4675270958133519, 1.4113417349322352, 3.44466700477302].+deltz)
]

3.9061001189111266e-7
x0 = [
    nSig(1, 1, 1.0349441753289228, 7.98002101545449e-5, Float64[], [0.0012497561308151344, 0.3462220260675253, 0.6207604655315174, 0.03176775227014218], [2.970265567210893, 4.834050801533279, 5.294620707583367, 6.5442719386535835])
    nSig(2, 1, 1.8037021012907866, 0.22288413646156122, [2.9493154300353353], [0.18563291134757007, 0.36030771750178486, 0.44433720945611943, 0.009722161694525673], [2.5211375527698063, 3.1437260844183683, 4.034094123935242, 6.098397008619959])
    nSig(3, 1, 1.9798994841803257, 0.19502503321097842, [0.5501562201339536, 2.9426027532622356], [0.006457064702259731, 0.0868737732502076, 0.4628588805958793, 0.4438102814516534], [0.23038448872030814, 1.067007067741543, 1.978260375750222, 3.6092504949765214])
    nSig(2, 2, 1.1778512161671362, 0.13698442149147613, Float64[], [0.0003391797880689138, 0.5118894702858586, 0.43306735910991806, 0.05470399081615445], [0.5647899716921982, 3.0867793180764433, 4.207843811259324, 5.825731425437815])
    nSig(3, 2, 1.893932026226742, 0.1913976248014796, [1.2521983341095697], [0.012954985369205153, 0.12633477793374373, 0.6246795349412988, 0.23603070175575236], [-0.15311136514912355, 0.6684568916990916, 1.5716470502433018, 3.5346789538327585])
]

Elim = -526.2745343;
ζ::Int16 = 18
ne::Int16 = 17
qt = [
    ETabn(-5,20,50),
    ETabn(-5,15,40)
]

x0 = [
    nSig(1, 1, 1.0367695364495866, 1e-3, Float64[], [0.001276837951486954, 0.34201177415856276, 0.6226922311756361, 0.03401915671431427], [2.7996198519362414, 4.7160286475148325, 5.182955163343529, 6.43005828877047])
    nSig(2, 1, 1.8164287483582775, 0.22280876454319756, [2.8344173833300093], [0.173196526504902, 0.3574556409679281, 0.4590691022249334, 0.01027873030223656], [2.362077836615896, 3.006583583628401, 3.903486390972554, 5.967283236044136])
    nSig(3, 1, 1.9277084437825913, 0.19501595158806834, [0.41324396270039504, 2.820390624753934], [0.0030266717530555588, 0.0533212494296452, 0.5176593839406324, 0.4259926948766669], [-0.123195973688624, 0.7034185946618937, 1.6733027359323076, 3.4262960763185886])
    nSig(2, 2, 1.175173203207956, 0.13800245967980745, Float64[], [8.934255954285922e-5, 0.48884977099607013, 0.4513331974614988, 0.05972768898288823], [0.16953043820732777, 2.9351988058236023, 4.063902885873775, 5.681658738455484])
    nSig(3, 2, 1.787623633012238, 0.21298255030688235, [1.0785040392562348], [0.008185227511168926, 0.10500337321483771, 0.6215390746633012, 0.2652723246106921], [-0.45769323473291124, 0.4675270958133519, 1.4113417349322352, 3.44466700477302])
]

3.6672668102255557e-7
x0 = [
    nSig(1, 1, 1.0370945326111591, 0.0009884543917008414, Float64[], [0.00123709476354829, 0.34068263732638776, 0.6240529280786936, 0.03402733983137044], [2.809675838706203, 4.717309694240429, 5.181698029112452, 6.431676794067673])
    nSig(2, 1, 1.8316887912487516, 0.22283627418388027, [2.8343591820719864], [0.1683622729270225, 0.3655133292858413, 0.4556900065396309, 0.010434391247505306], [2.3447433420432273, 2.9995004825533953, 3.9079845440927636, 5.963599607565825])
    nSig(3, 1, 2.102784411027585, 0.19012411671560342, [0.42187472458294134, 2.8165214981213826], [0.004881707526726137, 0.06146309939344473, 0.5644474320366341, 0.3692077610431951], [-0.118112735163727, 0.6185496903861725, 1.5831556902444288, 3.360070743923579])
    nSig(2, 2, 1.1515131845215676, 0.1278886273725246, Float64[], [5.589595210387803e-5, 0.49389755297344434, 0.44983312075838555, 0.05621343031606614], [0.41511540142748776, 2.9501934633539455, 4.073626181021137, 5.6966002062546535])
    nSig(3, 2, 2.0301716838493844, 0.21294584401119596, [1.084550959697847], [0.01430978990692124, 0.13438663907918622, 0.6043457037233236, 0.24695786729056907], [-0.3667152422758031, 0.455904188391276, 1.3538823349653102, 3.354622009680486])
]

Elim = -459.0485907;
ζ::Int16 = 17
ne::Int16 = 16
deltz = +2*log(17/18)
qt = [
    ETabn(-5+deltz,20,50),
    ETabn(-5+deltz,15,40)
]

x0 = [
    nSig(1, 1, 1.0379760188011171, 0.0007745008183693353, Float64[], [0.0013277417860542468, 0.32896470522843474, 0.6329315621823033, 0.03677599080320784], [2.5944750006915362, 4.590162356677897, 5.059580700545347, 6.301844780448613])
    nSig(2, 1, 1.832573071850657, 0.2227951052486209, [2.7123255323022883], [0.16405974525276534, 0.3612334362234311, 0.4636525674614534, 0.01105425106235015], [2.1941698192815706, 2.8783218116023135, 3.777375279286532, 5.8381510401107795])
    nSig(3, 1, 1.930038728414895, 0.19357284081068235, [0.23385625886345804, 2.7002808245991363], [0.0030961641534708565, 0.054059670397891, 0.489199272303119, 0.4536448931455191], [-0.30754225148209197, 0.5264582929676811, 1.4965633688061732, 3.271210996211579])
    nSig(2, 2, 1.1444577613925027, 0.1283197834183237, Float64[], [3.6847415659189036e-5, 0.4698860315095083, 0.4690748630847534, 0.061002257990079196], [0.04705211468342179, 2.7913266872729032, 3.924017199289512, 5.550129584514678])
    nSig(3, 2, 1.7145900848083555, 0.19204872578372426, [0.8891583524365174], [0.008601763382873408, 0.10279750049216747, 0.5966471684114231, 0.2919535677135361], [-0.6237740039802383, 0.2804519577094909, 1.2209354767667056, 3.315732923687774])
]

3.08004928228911e-7
x0 = [
    nSig(1, 1, 1.0379816747688586, 0.0007648398236103851, Float64[], [0.0013214383069096575, 0.3289709503566915, 0.632920070939586, 0.036787540396812826], [2.5945654565352507, 4.5903006099466195, 5.059466605287558, 6.301983717168589])
    nSig(2, 1, 1.8545716058427615, 0.22279348914269992, [2.7130684623666013], [0.1543682421370408, 0.36385419301152216, 0.4707980648260662, 0.010979500025370934], [2.1654272480934447, 2.849097856523261, 3.7659737716494655, 5.805483567380008])
    nSig(3, 1, 2.085584197547435, 0.17652710786094733, [0.24039762880881324, 2.6995869361428384], [0.004661894841141405, 0.06109764952602288, 0.5187636783026501, 0.4154767773301857], [-0.3040821241254864, 0.44376193683103626, 1.4295569225605727, 3.2158333882080026])
    nSig(2, 2, 1.133266198146665, 0.12371657137797552, Float64[], [3.5929288912458045e-5, 0.472211253434739, 0.4685029076457056, 0.05924990963064288], [0.27132145156842974, 2.7990507016973765, 3.9285653780560423, 5.556220838745285])
    nSig(3, 2, 2.1274705289716156, 0.22783395342132265, [0.8955029742094981], [0.01747146605242726, 0.14984390798705563, 0.5696629227115267, 0.2630217032489904], [-0.5495184003214263, 0.26325258777283345, 1.141306907206317, 3.182468383164567])
]

Elim = -397.1731828;
ζ::Int16 = 16
ne::Int16 = 15
deltz = +2*log(16/18)
qt = [
    ETabn(-5+deltz,20,50),
    ETabn(-5+deltz,15,40)
]

x0 = [
    nSig(1, 1, 1.0389367974960781, 3.225017958222763e-5, Float64[], [0.0010867587813183253, 0.3018408043333922, 0.6563358618779366, 0.04073657500735284], [2.4825148505955497, 4.441174168707965, 4.922814313949181, 6.1567896852870625])
    nSig(2, 1, 1.8484742734990611, 0.22285274979116207, [2.5831080511818554], [0.15185644005828103, 0.3605698297111777, 0.4760283640644101, 0.011545366166131187], [2.0108982972939993, 2.7287030229920486, 3.6304430441889206, 5.6838952344610245])
    nSig(3, 1, 1.962259038916618, 0.2033463597411732, [0.027611312081883853, 2.5752512456024483], [0.0035654478972354724, 0.06047420816039615, 0.439468623509643, 0.49649172043272544], [-0.48946974974128127, 0.3729040648921258, 1.3337468313985348, 3.1105764337008885])
    nSig(2, 2, 1.1206931740012231, 0.12282297625641332, Float64[], [2.5143332559903545e-5, 0.44645014491014123, 0.48961542765349614, 0.0639092841038027], [-0.07752728706895738, 2.6297791103846295, 3.7708040078256384, 5.403547091780809])
    nSig(3, 2, 1.7571127550129775, 0.19165012994325828, [0.6748441986755639], [0.009516169927671026, 0.10631018753446238, 0.5713109657980718, 0.3128626767397948], [-0.8655060434963158, 0.01750595507011512, 0.954046517892084, 3.1183573147395687])
]

3.1012564249977004e-7
x0 = [
    nSig(1, 1, 1.0389467020822192, 3.069369908635429e-12, Float64[], [0.0010641584187436157, 0.30202833234888776, 0.6562412102348403, 0.04066629899752843], [2.4915947572401325, 4.441160127367527, 4.922934637381182, 6.158037669665702])
    nSig(2, 1, 1.8626818007006567, 0.2225988910758943, [2.5833609398521555], [0.14724685036343252, 0.3639777045970291, 0.4772013554801688, 0.011574089559369644], [1.9941273159934838, 2.7152882460161436, 3.6285813299505607, 5.673427450069555])
    nSig(3, 1, 2.1492733688843697, 0.184864646895828, [0.03049378489459554, 2.5749596084383453], [0.005436648106192845, 0.06817187605196771, 0.4712463962036288, 0.4551450796382106], [-0.49944954625428817, 0.26878239352225997, 1.2628079768050757, 3.0477981479182117])
    nSig(2, 2, 1.1207789227590623, 0.12250451873933335, Float64[], [1.0272428489078879e-5, 0.44540867755529673, 0.4901145869320837, 0.06446646308413045], [0.0004925651483699481, 2.62841510045296, 3.767915603240599, 5.395664629154747])
    nSig(3, 2, 2.0337316515921184, 0.20603572864782754, [0.6720495548897235], [0.025679284631665397, 0.1776184315995336, 0.5008613100116953, 0.2958409737571058], [-0.6350097265180754, 0.1283279582111385, 0.9465908414934482, 3.024040950891729])
]

Elim = -340.3497759;
ζ::Int16 = 15
ne::Int16 = 14
deltz = +2*log(15/18)
qt = [
    ETabn(-5+deltz,20,50),
    ETabn(-5+deltz,15,40)
]

x0 = [
    nSig(1, 1, 1.03987891163068, 5.854474264301482e-5, Float64[], [0.000917505786830221, 0.29121929233756755, 0.6638825616845586, 0.04398064019104375], [2.519949069080526, 4.297030534336473, 4.785262776089519, 6.014444632244043])
    nSig(2, 1, 1.8487203794677864, 0.2202653600907479, [2.445409913075193], [0.13813506536697281, 0.353619987007336, 0.49596860634708134, 0.012276341278609914], [1.8190934771483385, 2.563880268846861, 3.46824246004038, 5.529661680533533])
    nSig(3, 1, 2.012863887961246, 0.21289193530639322, [-0.20173489943016948, 2.436626828124026], [0.004268422375392145, 0.06732065853434811, 0.417760723060313, 0.5106501960299468], [-0.7268866828510042, 0.15355356592831004, 1.1022144773909281, 2.909274028586538])
    nSig(2, 2, 1.0952747245423538, 0.11800929749281883, Float64[], [1.3805232877365579e-5, 0.4195923236205712, 0.5130553161135466, 0.06733855503300483], [-0.09231705338914949, 2.4527444929152495, 3.6060667789678686, 5.247067870903506])
    nSig(3, 2, 1.8211387323094952, 0.19169911319773472, [0.42044314644523634], [0.012544176696371687, 0.12302917765737101, 0.5255055541947613, 0.3389210914514961], [-1.0961833724332137, -0.2539505236423056, 0.6575568349044694, 2.902101241972619])
]

3.397261139070906e-7
x0 = [
    nSig(1, 1, 1.0399764421934614, 1.7573371007659022e-5, Float64[], [0.000918591732310043, 0.2912862819462764, 0.6637056220639426, 0.044089504257470984], [2.5729566182277477, 4.297191333219085, 4.784992725995354, 6.014116596254624])
    nSig(2, 1, 1.8629205662996924, 0.21980815914993648, [2.4452459271427744], [0.13621618438933236, 0.3615512692075821, 0.490010182055061, 0.01222236434802452], [1.805834178558763, 2.5617478928445316, 3.4748789663577813, 5.531798404215711])
    nSig(3, 1, 2.2984085459022454, 0.21280302950504776, [-0.2115129113634367, 2.4387921180820307], [0.006999008463080486, 0.08590764829764837, 0.4323247035079221, 0.47476863973134903], [-0.733396115014275, 0.06585350513716878, 1.0610106585143697, 2.838703835787371])
    nSig(2, 2, 1.0955436090656832, 0.11775631063021577, Float64[], [1.7424026063985244e-5, 0.41842042972224774, 0.5135492886003844, 0.06801285765130391], [-0.15263138624237985, 2.4512653073995736, 3.6028071673467856, 5.238697929246916])
    nSig(3, 2, 1.9979470225057903, 0.1916908637728561, [0.43280996237151825], [0.026193093992747864, 0.15993878281524826, 0.4969853274892859, 0.316882795702718], [-0.8781753247226599, -0.1694513633893486, 0.6545236333314689, 2.8564402766715937])
]

Elim = -288.5731311;
ζ::Int16 = 14
ne::Int16 = 13
deltz = +2*log(14/18)
qt = [
    ETabn(-5+deltz,20,50),
    ETabn(-5+deltz,15,40)
]

x0 = [
    nSig(1, 1, 1.0421129372312892, 3.63067808850405e-5, Float64[], [0.0009335212915466466, 0.2878102671375756, 0.6634704064104575, 0.047785805160420164], [2.4707614365327513, 4.147683942527546, 4.6410071917148965, 5.867033330959093])
    nSig(2, 1, 1.7416215799115, 0.18829226355721726, [2.298437287566229], [0.11496366235210545, 0.3208744828718413, 0.5517099630918801, 0.01245189168417324], [1.6321382250844427, 2.3566289008279155, 3.269117563910751, 5.428724323107054])
    nSig(3, 1, 1.7518644540259203, 0.15880240395599354, [-0.45355536674277414, 2.2826584931787197], [0.0037667960119768945, 0.05011864580623497, 0.3921287591600975, 0.5539857990216907], [-0.9088934777772883, -0.09439070733563477, 0.8614982205963848, 2.7460460901973964])
    nSig(2, 2, 1.0644805392890795, 0.11293126949098097, Float64[], [1.6140098082610655e-5, 0.38917690273831207, 0.5400910966578224, 0.07071586050578296], [-0.43697661715203195, 2.2596761818805904, 3.429605616188238, 5.082194143373481])
    nSig(3, 2, 1.7977971955654237, 0.19137576593581895, [0.14069517854253008], [0.014553624828901114, 0.12270445295936643, 0.4958062948963386, 0.3669356273153939], [-1.3826977768813833, -0.5413072405190582, 0.37090817925929526, 2.7456986957612948])
]

4.4168632484797854e-7
x0 = [
    nSig(1, 1, 1.042066093226808, 8.90653819341568e-6, Float64[], [0.0009355193319497013, 0.28748547145288883, 0.663758475931553, 0.047820533283608566], [2.492955948300853, 4.147784774561737, 4.640676084419606, 5.866817913271235])
    nSig(2, 1, 1.6515350624170528, 0.165614784401995, [2.296959026839758], [0.11722794365406225, 0.3164265634211069, 0.5532192167054121, 0.013126276219418647], [1.6776773167480217, 2.3864520285819038, 3.2864864739653057, 5.518399952431658])
    nSig(3, 1, 2.015996527026135, 0.17251176723095402, [-0.4581744824954632, 2.279609261372786], [0.010318833300096317, 0.08382809981096379, 0.42120736063740566, 0.4846457062515343], [-0.8212168683457093, -0.08923959571831488, 0.8042954916258697, 2.645831470361398])
    nSig(2, 2, 1.0640329896225793, 0.11260498844734197, Float64[], [1.1956765750979392e-5, 0.3883144045622376, 0.5405173106353679, 0.07115632803664354], [-0.4226586681950322, 2.259081973620751, 3.4272206724660466, 5.075915266778182])
    nSig(3, 2, 2.0898262908512146, 0.17690723700742056, [0.15598447290652567], [0.023809564103274214, 0.16345552049135537, 0.483724995244233, 0.3290099201611375], [-1.2695402867157894, -0.5472122875007002, 0.32855033605037925, 2.6602566217774233])
]

Elim = -241.6746705;
ζ::Int16 = 13
ne::Int16 = 12
deltz = +2*log(13/12)
qt = [
    ETabn(-5.2+deltz,20,50),
    ETabn(-5.2+deltz,15,40)
]

x0 = [
    nSig(1, 1, 1.0441225624046762, 3.533905549659831e-5, Float64[], [0.0008334600904358034, 0.2711262264499317, 0.6752435892167834, 0.052796724242849105], [2.252925362179744, 3.976786277518209, 4.479247564626897, 5.702035652920765])
    nSig(2, 1, 1.7200069887709557, 0.17426745308778693, [2.1395014253467406], [0.09106280111923455, 0.30601294289353886, 0.5901149939384547, 0.012809262048772028], [1.3739112110201441, 2.1253174085820903, 3.067789129436835, 5.255463695347774])
    nSig(3, 1, 1.68426341916748, 0.15853783294144289, [-0.7832198116731512, 2.126439328455716], [0.0046987213980178456, 0.0593464450068188, 0.30307180946322954, 0.6328830241319339], [-1.1188335377161187, -0.2415621107809719, 0.7130044139975036, 2.5763346836853938])
    nSig(2, 2, 1.0279512992270843, 0.10837288775756515, Float64[], [1.196033889283147e-5, 0.35417423941276477, 0.5719562868115664, 0.07385751343677596], [-0.7327980428974099, 2.046333508018852, 3.239850582218572, 4.908559405431543])
    nSig(3, 2, 1.853595434637422, 0.19087411735702858, [-0.20740190949386114], [0.014294619559342246, 0.11356379549541494, 0.4857142848636275, 0.3864273000816153], [-1.8691925406031544, -1.0163238050682644, -0.044829965406499875, 2.524846408787705])
]

5.569062011545611e-7
x0 = [
    nSig(1, 1, 1.04364606198283, 3.381114781171583e-5, Float64[], [0.000882760026831773, 0.2729728265262845, 0.6734746700783181, 0.052669743368565655], [2.3227321861646923, 3.980107487017373, 4.479855608618336, 5.701076364186155])
    nSig(2, 1, 1.5405629099767542, 0.13509177129376593, [2.134946297481322], [0.0991758954161299, 0.2995844965966348, 0.5855640995556045, 0.015675508431630777], [1.488549752102773, 2.2120619926038994, 3.1207947615566236, 5.417001732492266])
    nSig(3, 1, 1.9710419214146921, 0.15800251209156735, [-0.7807574493488364, 2.1208137687727175], [0.010741326494939075, 0.07786558817151758, 0.3699702209342724, 0.541422864399271], [-1.106695915119043, -0.4064895107517983, 0.5654699651353178, 2.427947564322257])
    nSig(2, 2, 1.021852486868064, 0.10651395959037922, Float64[], [3.5783764261422153e-6, 0.3549828162076576, 0.5726302532473407, 0.07238335216857562], [-0.8034734760755108, 2.0522537854860885, 3.2428946324566064, 4.9150035541550565])
    nSig(3, 2, 2.1439938791301842, 0.19087411735702858, [-0.23956393247278318], [0.017754198979864258, 0.1131721856461558, 0.5095495670775871, 0.35952404829639284], [-2.1532560578612565, -0.9978519771611598, -0.04461721111938592, 2.6322003845509947])
]

Elim = -199.3718097;
ζ::Int16 = 12
ne::Int16 = 11
qt = [
    ETabn(-5.2,20,50),
    ETabn(-5.2,15,40)
]

x0 = [
    nSig(1, 1, 1.0462809137163462, 1.3315767901114318e-11, Float64[], [0.0010428899424284306, 0.28417367549439865, 0.6579570325428103, 0.05682640202036258], [2.24515992878133, 3.819418728389465, 4.318683355116902, 5.535983190151826])
    nSig(2, 1, 1.651587064594959, 0.16231760430062694, [1.9640832945076054], [0.07930776805849307, 0.29508618353400384, 0.611099214106165, 0.014506834301338244], [1.1529732220582074, 1.9452335191599162, 2.8852036114703083, 5.146880149311541])
    nSig(3, 1, 1.5812018918236375, 0.12792925036360947, [-1.1811561075504513, 1.9438093801589], [0.005125853244088781, 0.0477967746976934, 0.31096247110743597, 0.6361149009507819], [-1.548569519918797, -0.7006115511560769, 0.3131987582033363, 2.2973041124840874])
    nSig(2, 2, 0.985055199807286, 0.10532806224967234, Float64[], [0.14, 0.36, 0.40, 0.10], [1.2, 2.2, 3.2, 4.7])
]

4.4845342017652e-7
x0 = [
    nSig(1, 1, 1.0462669849601667, 9.560165969626827e-13, Float64[], [0.0010490899339422248, 0.28415531462821764, 0.6579754550857111, 0.05682014035212911], [2.2479202490095367, 3.8193855022816545, 4.318722421958952, 5.535998467860176])
    nSig(2, 1, 1.9299451098649185, 0.24561648147045867, [1.9638612021302866], [0.0923613043544756, 0.3496491484480367, 0.5428648190955895, 0.015124728101898278], [1.0960767946008547, 1.982656256220365, 2.924608922472219, 4.996007033952849])
    nSig(3, 1, 1.2741619559474997, 0.07466739567727282, [-1.1806160705569477, 1.939864134805536], [0.029890964095886513, 0.04119501073320901, 0.28591304598058337, 0.6430009791903212], [-0.6164997334539084, -0.601512177879806, 0.29387558150173415, 2.3977410479450443])
    nSig(2, 2, 1.1416582941874296, 0.11863264471910802, Float64[], [0.19134204484430525, 0.3491305659545624, 0.37255962146490385, 0.08696776773622855], [1.456561467839866, 2.4106320835415076, 3.2322643743375594, 4.77552372038184])
]

Elim = -161.6769626;
ζ::Int16 = 11
ne::Int16 = 10
deltz = +2*log(11/10)
qt = [
    ETabn(-4.2+deltz,20,50),
    ETabn(-4.2+deltz,15,40)
]

x0 = [
    nSig(1, 1, 1.0500762892445639, 0.0007801787179705028, Float64[], [0.0010775167907148683, 0.29533809769109653, 0.64202261197174, 0.061561773546448574], [2.0588832569206392, 3.6430923128528505, 4.142567206060774, 5.360115268314258])
    nSig(2, 1, 1.5574736157611417, 0.16235678387802108, [1.7693974496510043], [0.07547365601663444, 0.29712687224432316, 0.6076757081702294, 0.019723763568813085], [0.9692906205669075, 1.8267318806846522, 2.729481657241819, 5.033097902067228])
    nSig(2, 2, 1.197430214241464, 0.15014520042507845, Float64[], [0.1404418545419179, 0.36158349021291564, 0.3948635648946892, 0.10311109035047726], [1.0309070879349842, 2.0797635720391505, 2.996479179271747, 4.55826651574818])
]

1.2746983202305273e-7
x0 = [
    nSig(1, 1, 1.0497508180202197, 3.177418796198292e-10, Float64[], [0.0011356959677671032, 0.2964999131031466, 0.640584288030913, 0.061780102898173234], [2.084481406470743, 3.643038164802748, 4.14302346244602, 5.3571921405185705])
    nSig(2, 1, 1.5622217740400395, 0.1408289133332779, [1.771804154787943], [0.06361742131573886, 0.27933413100042137, 0.6403323992520912, 0.016716048431748617], [0.88776991497044, 1.7341484257583253, 2.679777357442932, 5.005175724216771])
    nSig(2, 2, 1.0728553225199662, 0.06744140602944765, Float64[], [0.12112121545921697, 0.3241427172846334, 0.4573419234816042, 0.09739414377454549], [0.9947905763303126, 1.9745517267169017, 2.9329353475419957, 4.520182474982323])
]

Elim = -127.8178141;
ζ::Int16 = 10
ne::Int16 = 9
qt = [
    ETabn(-4.2,20,50),
    ETabn(-4.2,15,40)
]

x0 = [
    nSig(1, 1, 1.0526107361511146, 0.0007731119007157296, Float64[], [0.0012536240619454746, 0.309171850553032, 0.6219872778122594, 0.06758724757276326], [1.9397025364531264, 3.4538942865137416, 3.948384184664512, 5.157950737033293])
    nSig(2, 1, 1.5053371219704514, 0.1545608103406165, [1.5598628380555148], [0.04707225978041301, 0.2592149083081062, 0.6721034614351898, 0.02160937047629098], [0.6017566482148391, 1.4975653094516208, 2.459822365409779, 4.823441547156512])
    nSig(2, 2, 1.0397866768729789, 0.08601803577726749, Float64[], [0.07805723864227287, 0.2982461789169043, 0.5137815241101246, 0.10991505833069819], [0.5295334914961968, 1.6125709313049617, 2.64643804044916, 4.2844190203535675])
]

7.844728600048256e-8
x0 = [
    nSig(1, 1, 1.0531235840510786, 0.00076736571047296, Float64[], [0.0011099568602472708, 0.3064358081643447, 0.6245299405551136, 0.06792429442029445], [1.9074534140752926, 3.4516782085581283, 3.9460457770345614, 5.159019297204496])
    nSig(2, 1, 1.6048411848697222, 0.1436864099994203, [1.5639934355926244], [0.06593697878911514, 0.2801793892363377, 0.6354924669029207, 0.018391165071626525], [0.6668874435132095, 1.502336620860351, 2.449201115825254, 4.789641651633049])
    nSig(2, 2, 1.0727637586002023, 0.060925889306696304, Float64[], [0.13262938038131092, 0.32263281620059564, 0.43788448287889514, 0.10685332053919822], [0.8301205803604736, 1.7795841339028202, 2.7230024679278335, 4.322181610354807])
]

Elim = -98.83172020;
ζ::Int16 = 9
ne::Int16 = 8
deltz = +2*log(9/10)
qt = [
    ETabn(-4.2+deltz,20,50),
    ETabn(-4.2+deltz,15,40)
]

x0 = [
    nSig(1, 1, 1.0550830588507665, 0.000754806724223887, Float64[], [0.0009780428485369236, 0.29413322946421405, 0.6284488881733814, 0.07643983951386762], [1.640143596453736, 3.2206951974886113, 3.7165994631056245, 4.926959574949433])
    nSig(2, 1, 1.4808627589866021, 0.14343079955930033, [1.3250067233354252], [0.04858655601411737, 0.25530466645676464, 0.6695450508864501, 0.02656372664266796], [0.3847709539282141, 1.2658559279290995, 2.2230834509242188, 4.606732720508366])
    nSig(2, 2, 1.0559811248517048, 0.08601531159010287, Float64[], [0.08492889390435907, 0.30201307155636115, 0.4883403591232095, 0.12471767541607026], [0.31213450776515267, 1.380426080038497, 2.403126970724038, 4.060102022688833])
]

6.076642478092253e-8
x0 = [
    nSig(1, 1, 1.0550309118939811, 0.0007532805664198324, Float64[], [0.00081981989912891, 0.28684829857542166, 0.6353236604629651, 0.07700822106248441], [1.5818075940581986, 3.2163902141346705, 3.711640406520453, 4.925390592454211])
    nSig(2, 1, 1.5733460100317018, 0.12915568040554062, [1.3292430898447047], [0.07034939216538805, 0.2753484560034535, 0.6319923788591978, 0.022309772971960677], [0.4683228624138012, 1.2731103672325634, 2.210733392031479, 4.594413886434473])
    nSig(2, 2, 1.0955033691948763, 0.061945014444065076, Float64[], [0.15268314487911808, 0.32982768969419146, 0.3969965453040123, 0.12049262012267814], [0.6545492362785684, 1.5801849405150363, 2.5041920855511592, 4.099057368607323])
]

Elim = -74.37260568;
ζ::Int16 = 8
ne::Int16 = 7
deltz = +2*log(8/10)
qt = [
    ETabn(-4.2+deltz,20,50),
    ETabn(-4.2+deltz,15,40)
]

x0 = [
    nSig(1, 1, 1.058506475477142, 0.0007553886688798232, Float64[], [0.000824259018831871, 0.2944697403091797, 0.6186834491536648, 0.08602255151832365], [1.346506967760071, 2.970603837669832, 3.465283402490668, 4.676868453931501])
    nSig(2, 1, 1.5083664541602089, 0.14341477964129953, [1.059592185427134], [0.05071392027394572, 0.25416153781318623, 0.6651954089279957, 0.029929132984872295], [0.10413878307782082, 0.978799254774404, 1.9342368646239674, 4.352625388528668])
    nSig(2, 2, 1.060138143081228, 0.07930945766734267, Float64[], [0.0970129650204456, 0.30673133440896305, 0.4546704907267403, 0.14158520984385115], [0.09761127957978243, 1.1289791047316098, 2.1329875261343787, 3.805636994087329])
]

4.316802915127482e-8
x0 = [
    nSig(1, 1, 1.0592085795518849, 8.468758516264851e-9, Float64[], [0.0006982495085670005, 0.293945155116634, 0.6191466084425444, 0.08620998693225453], [1.2928745598310944, 2.968333092734463, 3.464092778612546, 4.679826804588718])
    nSig(2, 1, 1.6694094743857029, 0.1380345331365614, [1.0659014463407661], [0.07182385217991984, 0.2748615550258361, 0.6305109701494346, 0.022803622644809517], [0.15891055377647306, 0.947496055403221, 1.8997158700047934, 4.321403025342009])
    nSig(2, 2, 1.1086955533927538, 0.06066875420960255, Float64[], [0.19463298647319177, 0.346792984390206, 0.32285158378216205, 0.13572244535444014], [0.5139458297347924, 1.4033441915290483, 2.3030385373965276, 3.8482929090219344])
]

Elim = -53.88800501;
ζ::Int16 = 7
ne::Int16 = 6
deltz = +2*log(7/10)
qt = [
    ETabn(-4.2+deltz,20,50),
    ETabn(-4.2+deltz,15,40)
]

x0 = [
    nSig(1, 1, 1.0634207753378446, 0.0007256102997371489, Float64[], [0.0008569030992198541, 0.32305532625158934, 0.5801053463628962, 0.0959824242862947], [1.1143253041281864, 2.703374392984457, 3.1948432604761225, 4.403299023477875])
    nSig(2, 1, 1.5497834857970576, 0.14288068539876023, [0.7519355669974236], [0.05290272045208332, 0.2499305917189128, 0.6640361443488666, 0.03313054348013729], [-0.2346326151413587, 0.6318046239614603, 1.592645550289613, 4.065534200548969])
    nSig(2, 2, 1.0704354770242508, 0.07152173189245406, Float64[], [0.1232691455467529, 0.3220148822557212, 0.39017931657817323, 0.16453665561935274], [-0.10008871291064717, 0.8665602572701611, 1.8456869346236926, 3.506294454228617])
]

2.6543467868123116e-8
x0 = [
    nSig(1, 1, 1.0631490481006665, 0.0007205241788514443, Float64[], [0.0005425118469885946, 0.3075871791291017, 0.594355029815912, 0.09751527920799771], [1.0203239165623603, 2.6964408059550697, 3.182659196416096, 4.39873743450065])
    nSig(2, 1, 1.7257888059091213, 0.1385941581255669, [0.7611315266844105], [0.07666795862021165, 0.27859930643415054, 0.619050727334509, 0.025682007611128873], [-0.1438801217301894, 0.6152032024606329, 1.5616604077818466, 4.04776556002543])
    nSig(2, 2, 1.1551047691715193, 0.06172626822130668, Float64[], [0.2089903932214732, 0.3424396243750886, 0.29552936710371086, 0.15304061529972737], [0.20842048433980215, 1.081609919072317, 2.004597211971281, 3.563847711668656])
]

Elim = -37.29222377;
ζ::Int16 = 6
ne::Int16 = 5
deltz = +2*log(6/10)
qt = [
    ETabn(-4.2+deltz,20,50),
    ETabn(-4.2+deltz,15,40)
]

x0 = [
    nSig(1, 1, 1.0668057585792166, 0.0007297142563670196, Float64[], [0.0007062421717345965, 0.33741295440475333, 0.5513967828837122, 0.11048402053979987], [0.8281511963901425, 2.3822536814848485, 2.868694678198442, 4.07449773154369])
    nSig(2, 1, 1.5906323264323772, 0.14286107544688062, [0.38735123854066233], [0.055588749858419315, 0.24608800317436869, 0.660562883611622, 0.037760363355590104], [-0.6126710666648524, 0.22717522381218985, 1.1846226071660273, 3.737684801968485])
    nSig(2, 2, 1.0696030039377376, 0.0663897427747004, Float64[], [0.14410664269445575, 0.31893597058650264, 0.34976513420373845, 0.1871922525153032], [-0.4337686026938916, 0.5237437620540063, 1.4840355186717593, 3.1816066302024857])
]

1.4294720074303768e-8
x0 = [
    nSig(1, 1, 1.0678674781398247, 0.00072779928987984, Float64[], [0.00032212036298679417, 0.3343407429170052, 0.5540494378968228, 0.11128769882318518], [0.6795876244991231, 2.382545282940387, 2.862368794252564, 4.076674238881732])
    nSig(2, 1, 1.7669163951078992, 0.13292752017795406, [0.39871820399414687], [0.08424143980017472, 0.28358249442563116, 0.6022436662063416, 0.02993239956785251], [-0.48132078367878955, 0.23250196824452868, 1.16340106571356, 3.7406246241668653])
    nSig(2, 2, 1.1249715715814923, 0.05191690358091854, Float64[], [0.27649484200940755, 0.32454250883455427, 0.22848681489789308, 0.17047583425814503], [-0.014786663123203314, 0.8568751325429051, 1.7154316032496837, 3.24805130464172])
]

Elim = -24.23757518;
ζ::Int16 = 5
ne::Int16 = 4
deltz = +2*log(5/4)
qt = [
    ETabn(-5.4+deltz,20,50)
]

x0 = [
    nSig(1, 1, 1.071190841673446, 0.0006814743068373455, Float64[], [0.0001, 0.35583874500921625, 0.514645401468512, 0.1291538997892242], [0.8, 1.9931404327420037, 2.483338969156536, 3.6879709111478634])
    nSig(2, 1, 1.1825128167043235, 0.04990260974322443, [-0.07233392468645437], [0.0711522145995949, 0.22398761776503112, 0.5983340880621957, 0.10652607957317832], [-0.6937180268187345, 0.06069545265126659, 0.9908995098213567, 3.253983878176012])
]
# Substitui o x0[2] mais uma vez e otimiza só ele 2x.

5.203258268693389e-9
x0 = [
    nSig(1, 1, 1.0706612049269746, 1.236080759871254e-6, Float64[], [0.0008043522247147328, 0.36022802286644473, 0.5083587983901823, 0.13060882651865838], [1.8259198389719746, 2.00413490528556, 2.476702496420774, 3.686180966966906])
    nSig(2, 1, 1.1633606425135625, 0.014825070853866083, [-0.058213315373148555], [0.12210127882308133, 0.25710229511360566, 0.5125282356944563, 0.10826819036885667], [-0.3705406827825325, 0.15190290498645512, 1.0743770771049694, 3.267620738089765])
]

Elim = -14.27739481;
ζ::Int16 = 4
ne::Int16 = 3
qt = [
    ETabn(-5.4,20,50),
]

x0 = [
    nSig(1, 1, 1.0755274242531838, 0.0006805001470837334, Float64[], [0.0005242406872370283, 0.4136154715612843, 0.4339804248554629, 0.1518798628960157], [0.5783125646449293, 1.5251683145457324, 2.0406066091328747, 3.219815717815142])
    nSig(2, 1, 1.1574154646988064, 0.03536521208859407, [-0.6685171882083422], [0.09241347893203239, 0.2104689461297181, 0.5600412515245033, 0.13707632341374618], [-1.2157162193175877, -0.5220919194473091, 0.40729936533205424, 2.750410189121941])
]

2.1966872765233347e-9
x0 = [
    nSig(1, 1, 1.074295682577798, 0.0006426018566320401, Float64[], [0.00020571090042634856, 0.456398476828308, 0.39322628734735104, 0.15016952492391467], [-0.8974840813524155, 1.5496703056469499, 2.07462177185157, 3.218270312732805])
    nSig(2, 1, 1.288818761042698, 0.03511252926123416, [-0.6633882690901995], [0.20572726517388965, 0.3255249669236132, 0.38201815864862193, 0.0867296092538752], [-0.8463896533794325, -0.5053367027132867, 0.12283188469987405, 2.7412787963091914])
]

Elim = -7.236415201;
ζ::Int16 = 3
ne::Int16 = 2
deltz = +2*log(3/2)
qt = [
    ETabn(-6+deltz,20,50),
]

x0 = [
    nSig(1, 1, 1.0773634222467896, 0.0006911243375054187, Float64[], [8.819644833281902e-5, 0.45828547227771965, 0.3613883889516055, 0.18023794232234197], [-1.9629358000331554, 0.8739436185611957, 1.493634293328408, 2.6165866859409763])
]

1.9591883670955212e-9
x0 = [
    nSig(1, 1, 1.0806501475839982, 0.0006911243375054187, Float64[], [0.0022069559776462925, 0.5354438651947679, 0.298843524736768, 0.16350565409081785], [0.19789328412983342, 0.9182376320761211, 1.6264048218353904, 2.6512970114268897])
]

Elim = -2;
ζ::Int16 = 2
ne::Int16 = 1
qt = [
    ETabn(-6+deltz,20,50),
]

x0 = [
    nSig(1, 1, 1, 0.01, Float64[], [1], [1.7])
]

1.9100276915651193e-12
x0 = [
    nSig(1, 1, 1.0000008587125695, 2.4840494449904403e-8, Float64[], [1.0], [0.6931465942723263])  
]

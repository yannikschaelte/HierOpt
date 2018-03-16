function [ D ] = getData(index)
    
    switch(index)  
        case 1
            
            t{1}  = [2 3 4 5 6 38 81 144 211 408 601 804 1204 2004]';
            t{2}  = [2 3 4 5 6 38 81 144 211 408 601 804 1204 2004]';
            t{3} = [2 3 4 5 6 38 81 144 211 408 601 804 1204 2004]';
            t{4}  = [2 3 4 5 6 38 81 144 211 408 601 804 1204 2004]';
            t{5} = [38 91 208 408 608 1008 1208 1418 1608 1804 2007]';
            t{6}  = [14 108 128 208 381 611 1004 1994]';
            t{7}  = [24 28 81 83 141 173 201 216 391 394 799 828 1188 1201 1394]';
            t{8}  = [31 74 134 201 401 601 798 1201 1598 2001]';
            t{9}  = [31 74 134 201 401 601 798 1201 1598 2001]';
            
            D.t = sort(unique([t{1};t{2};t{3};t{4};t{5};t{6};t{7};t{8}]));
            
            exp_data{1} = [ ...
                1.99;
                2.1;
                2.09;
                1.84;
                2.31;
                2.76;
                3.05;
                2.42;
                2.23;
                2.52;
                2.81;
                2.71;
                2.71;
                2.7;
                ];
            
            exp_data{2} = [ ...
                4.39;
                4.76;
                4.86;
                4.65;
                4.75;
                5.52;
                5.86;
                4.39;
                3.6;
                3.83;
                4.3;
                4.05;
                3.27;
                3.38;
                ];
            
            exp_data{3} = [ ...
                4.07;
                3.71;
                3.19;
                3.57;
                3.14;
                2.38;
                3.71;
                3.19;
                5.24;
                4.47;
                3.62;
                3.62;
                2.86;
                2.4;
                ];
            
            exp_data{4} = [ ...
                0.62;
                0.66;
                0.74;
                0.62;
                0.75;
                0.92;
                1.15;
                0.57;
                0.46;
                0.57;
                0.57;
                0.69;
                0.46;
                0.46;
                ];
            
            exp_data{5}	= [ ...
                1.255555556;
                1.311111111;
                1.283333333;
                0.8611111111;
                0.5972222223;
                0.09611111112;
                0.04333333334;
                0.05055555556;
                0.04777777778;
                0.04777777778;
                0.06;
                ];
            
            exp_data{6}	= [ ...
                1.35;
                0.83;
                0.83;
                0.78;
                0.84;
                0.64;
                0.74;
                0.70;
                ];
            
            exp_data{7}	= [ ...
                1.01;
                0.92;
                1.15;
                1.19;
                1.06;
                1.10;
                1.05;
                1.08;
                0.97;
                1.01;
                0.92;
                0.89;
                0.74;
                0.88;
                0.80;
                ];
            
            exp_data{8}	= [ ...
                0.19;
                0.56;
                1;
                2.83;
                1.5;
                2.26;
                2.4;
                1.25;
                0.07;
                0.02;
                ];
            
            exp_data{9}	= [ ...
                0.19;
                0.56;
                1;
                2.83;
                1.5;
                2.26;
                2.4;
                1.25;
                0.07;
                0.02;
                ];
            
            Y = NaN(length(D.t),9);
            
            for iobs=1:9
                for it = 1:length(t{iobs})
                    idx = find(t{iobs}(it)==D.t,1,'first');
                    Y(idx,iobs) = exp_data{iobs}(it);
                end
            end
            
            D.Y = Y;
            D.Sigma_Y = 0.15*ones(size(Y));
            
            D.condition = [];
            
        case 2
            
            load('b3_data.mat');
            D.Y = xnom;
            
            expDataMax = max(xnom);
            
            Q = repmat(1./expDataMax, size(D.Y,1), 1);
            atol = 1e-7;
            tmp1 = or(isnan(Q),Q > 1/atol);
            Q(tmp1) = 1;
            
            D.Sigma_Y = 1./Q;
            
            D.t = transpose([0:1000:159480,159480]);
            
            D.condition = [];
            
        case 3
            
            css(1)=100000;
            css(2)=1;
            css(3)=1000;
            css(4)=1000;
            css(5)=30;
            css(6)=1000;
            css(7)=1000;
            css(8)=1000;
            css(9)=1000;
            css(10)=1000;
            css(11)=1000;
            css(12)=1000;
            css(13)=3000;
            css(14)=1000;
            css(15)=1000;
            css(16)=1000;
            css(17)=100;
            css(18)=100;
            css(19)=100;
            css(20)=100;
            css(21)=1000;
            css(22)=1000;
            css(23)=1000;
            css(24)=1000;
            css(25)=1000;
            css(26)=1000;
            css(27)=1000;
            css(28)=1000;
            css(29)=1000;
            css(30)=3000;
            css(31)=100000;
            css(32)=1000;
            css(33)=1000;
            css(34)=1000;
            
            t_s      = [23.9759 48.4337 72.8916 90.3614 115.984 140.442 164.9 210.321 233.614 258.072 282.53 300];
            timeseries_artificalData_noised = [
                % Productprotein	L-Methionine	L-Leucine	L-Lactate	beta-D-Glucose	L-Aspartate	L-Malate	Pyruvate	Oxaloacetate	ATP	ATP	ADP	ADP_c
                % Time	f	f	f	f	f	c	c	c	m	c	m	m	c
                % min	mumol/l_fermenter	mumol/l_fermenter	mumol/l_fermenter	mumol/l_fermenter	mumol/l_fermenter	mumol/l_cell	mumol/l_cell	mumol/l_cell	mumol/l_cell	mumol/l_cell	mumol/l_cell	mumol/l_cell	mumol/l_cell
                0	    29.54174187	5110.119691	4999.956639	1.056682831	449958.4343	1009.945764	988.3517969	996.1405443	839.7703035	3314.558385	2968.296835	1016.447754	1029.663149
                23.9759	29.06135825	4427.537773	4990.471667	19033.05068	499628.9405	514.889607	1883.981851	4115.730274	443.5222307	3298.082924	3588.716676	701.1769714	841.6656476
                48.4337	31.56920926	4649.633023	4706.446115	29713.37459	498137.7445	614.0882509	1999.234422	5088.17481	457.2654501	3144.678363	3524.003949	534.3245626	898.4469128
                72.8916	32.58992629	4642.175831	4703.420191	52653.84782	469436.4399	667.144403	1977.174384	4941.17988	452.1008571	3343.981672	3274.933898	637.1604048	748.7508252
                90.3614	32.93492455	4581.482695	4580.779761	62491.60793	458816.0559	706.4697945	2027.556359	5561.084132	416.1288178	3305.264733	3377.224218	608.6810043	760.8237877
                115.984	33.80995822	4231.998041	4576.568758	71665.49964	441725.8779	748.837145	1834.640732	4711.566421	433.6235068	2664.13707	3359.281519	608.5031107	746.0614226
                140.442	34.67192676	4609.073999	4388.444795	93744.69241	371353.6441	761.4219172	1741.34438	4466.968205	471.9681186	3276.842196	3268.007202	609.2855454	770.8013591
                164.9	35.64086504	4379.397091	4823.258118	96012.12538	409872.6888	913.9581526	2054.553359	4841.62543	462.4384044	3326.232934	3375.331315	657.6144762	747.1622156
                210.321	35.1073701	4047.37735	4757.799632	131633.9417	429309.5792	838.4429003	2023.441102	4741.87247	463.3431076	2838.535673	3344.113884	606.0191858	696.1165299
                233.614	34.78766884	4354.532227	4267.505383	139701.0391	377830.4707	864.115653	2012.727053	4982.305375	451.6411544	3445.146289	3384.49464	596.595147	747.4997328
                258.072	35.10275329	3893.078648	3946.812817	156106.4324	376607.9654	1019.41341	2077.935601	4932.876021	432.2633547	3185.845081	3123.317817	606.9552722	742.2563724
                282.53	39.16229127	4223.934695	4229.786474	189350.8377	350136.0051	945.7524046	2018.02482	4965.805674	459.569613	3248.126308	3303.223152	634.2843344	740.5121927
                300	39.36873083	3819.952372	3226.591797	170253.9838	362466.0386	962.815767	2167.954632	4922.316021	461.1326466	3189.167366	3335.671338	664.509921	757.3732437];
            
            unscaled_data     = timeseries_artificalData_noised(:,2:end);
            scaled_data       = unscaled_data;
            scaled_data(:,1)  = unscaled_data(:,1)/css(5);
            scaled_data(:,2)  = unscaled_data(:,2)/css(4);
            scaled_data(:,3)  = unscaled_data(:,3)/css(3);
            scaled_data(:,4)  = unscaled_data(:,4)/css(2);
            scaled_data(:,5)  = unscaled_data(:,5)/css(1);
            scaled_data(:,6)  = unscaled_data(:,6)/css(29);
            scaled_data(:,7)  = unscaled_data(:,7)/css(27);
            scaled_data(:,8)  = unscaled_data(:,8)/css(21);
            scaled_data(:,9)  = unscaled_data(:,9)/css(15);
            scaled_data(:,10) = unscaled_data(:,10)/css(13);
            scaled_data(:,11) = unscaled_data(:,11)/css(30);
            scaled_data(:,12) = unscaled_data(:,12)/css(32);
            scaled_data(:,13) = unscaled_data(:,13)/css(11);
            
            timeseries_artificalData_stdev = [
                0	    0.916516259	220.2393825	0.086722941	0.113365662	100083.1314	19.8915281	23.29640617	7.718911471	320.4593929	629.11677	63.40633065	32.89550786	59.32629796
                23.9759	3.42708349	859.4844544	266.3833335	524.0986472	27089.88103	72.52278599	1.363702921	13.17945154	67.9675385	110.2258472	512.6933529	67.11194286	169.2792952
                48.4337	0.072781474	197.1939541	83.56777088	13530.25082	50101.48902	8.871498249	109.6688443	1402.38962	25.79309975	208.3432738	361.1878982	244.5448747	294.6078256
                72.8916	0.345652583	11.22833874	111.2603823	145.6956478	17694.8798	11.00519392	1.408768357	800.1597591	29.84428584	184.2633444	149.0122046	26.80319038	1.191650429
                90.3614	0.103150906	5.585389682	4.179521469	2222.584132	13862.11186	0.953588928	72.01271842	1881.128264	99.10036436	104.3694651	49.86843658	78.08799131	27.79957543
                115.984	0.003516434	495.8039179	193.3375168	15030.00072	4663.755834	0.759709907	343.8585363	3.332842875	61.44498639	1179.74586	8.403037386	72.83977854	0.150845275
                140.442	0.186053515	442.6879971	1.429589001	515.5848117	112764.7117	48.67616557	548.4312396	616.2835899	16.98823721	45.18439211	177.4455966	67.99290924	50.07871826
                164.9	0.604730081	164.3741824	1052.096236	22589.74924	12874.62245	190.2563052	67.32671764	29.79086077	0.791191271	144.5058684	35.34263088	30.52495245	2.300431195
                210.321	3.228859795	170.8252992	1250.019264	372.1165486	67335.15833	70.03419938	0.977796603	311.3950597	2.654215273	828.1886549	27.87223266	71.86362837	102.5249403
                233.614	5.262062329	608.8844539	434.830766	8311.921766	14941.05855	69.51269399	20.40589406	115.4107506	20.1576913	387.132577	53.48928055	91.2997059	1.854534475
                258.072	6.079093429	142.482705	35.01436682	49.13512714	3973.930766	190.684819	114.9912025	29.80795876	58.42129056	128.8298377	467.6643662	71.7674556	14.94525519
                282.53	0.608382544	688.7293909	700.4329474	42607.67548	27959.98971	4.487190854	2.669639537	1.328651269	3.420773959	1.32738423	106.1736951	18.79533114	21.39761464
                300	    0.008461654	0.624744428	1186.096406	12180.03246	11498.07728	3.204466046	309.209264	110.2279582	0.074706795	116.9052678	39.83732302	40.18584194	10.01248736];
            
            unscaled_error    = timeseries_artificalData_stdev(:,2:end);
            scaled_error      = unscaled_error;
            scaled_error(:,1)  = unscaled_error(:,1)/css(5);
            scaled_error(:,2)  = unscaled_error(:,2)/css(4);
            scaled_error(:,3)  = unscaled_error(:,3)/css(3);
            scaled_error(:,4)  = unscaled_error(:,4)/css(2);
            scaled_error(:,5)  = unscaled_error(:,5)/css(1);
            scaled_error(:,6)  = unscaled_error(:,6)/css(29);
            scaled_error(:,7)  = unscaled_error(:,7)/css(27);
            scaled_error(:,8)  = unscaled_error(:,8)/css(21);
            scaled_error(:,9)  = unscaled_error(:,9)/css(15);
            scaled_error(:,10) = unscaled_error(:,10)/css(13);
            scaled_error(:,11) = unscaled_error(:,11)/css(30);
            scaled_error(:,12) = unscaled_error(:,12)/css(32);
            scaled_error(:,13) = unscaled_error(:,13)/css(11);
            
            D.Y = scaled_data(2:end,:);
            D.Sigma_Y = scaled_error(2:end,:);
            
            D.t = t_s;
            
            D.condition = [];
            
        case 4
            
            D(1).t = [0  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30];
            
            D(1).Y = [...
                0.11327346 0.15762948 0.16568655 0.9004708 0.10767257 0.1355827
                0.19829967 0.14064558 0.17512615 0.9041830 0.09441348 0.2145996
                0.09299877 0.06633699 0.10856162 0.9375022 0.14480097 0.1825351
                0.12164888 0.11327575 0.15403017 0.9622286 0.10224031 0.1371931
                0.21959824 0.17236626 0.22072767 0.8919632 0.13063334 0.1498756
                0.15472042 0.12708298 0.11575208 0.8979801 0.22518701 0.1550437
                0.11013992 0.11169585 0.16922386 0.8078341 0.10325274 0.1670639
                0.18678235 0.05988388 0.12707947 0.8855724 0.11136211 0.1085992
                0.11089965 0.14960227 0.00000000 0.9211519 0.09265214 0.1506382
                0.13010459 0.05472079 0.14316334 0.8873765 0.15153122 0.1567842
                0.12903695 0.12898450 0.13989877 0.8459586 0.06496903 0.1988740
                0.12996955 0.18797274 0.11880129 0.8869094 0.24472758 0.1612959
                0.17110675 0.24895319 0.08181443 0.9426594 0.17549157 0.1795871
                0.07539216 0.08494266 0.13547178 0.9059816 0.22607564 0.1557416
                0.08757456 0.07813579 0.11815741 0.8782907 0.22484167 0.1848494
                0.09236308 0.14698945 0.15508290 0.8925849 0.22967481 0.1172746];
            
            D(1).Sigma_Y = 0.05*ones(16,6);
            
            D(1).condition = [0    0    0    0];
            
            D(2).t = [0  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30];
            
            D(2).Y = [...
                0.1612526 0.1361274 0.1321657 0.97346347 0.1260189 0.13759391
                0.6024987 0.7634884 0.3054694 0.64127482 0.1526805 0.18244736
                0.8575352 0.8334763 0.6127695 0.46434514 0.1657532 0.14044900
                0.7839275 0.8719813 0.5913194 0.25160181 0.2118012 0.16972103
                0.7449103 0.8196278 0.7289923 0.24798320 0.1734451 0.13354375
                0.7449549 0.8004981 0.8416676 0.25242939 0.1880558 0.09318322
                0.5856446 0.8208235 0.8224108 0.11910233 0.1511175 0.15746581
                0.6331345 0.7379835 0.8530437 0.15276035 0.1375105 0.12900267
                0.4925625 0.5589710 0.8763375 0.21345420 0.1274942 0.16265159
                0.4620047 0.4529328 0.7517932 0.12449379 0.1184203 0.15001744
                0.3258993 0.3197356 0.5493157 0.15283422 0.1177744 0.20602413
                0.2970832 0.2839845 0.3524049 0.08166205 0.1476658 0.14856343
                0.3005535 0.2521285 0.2718919 0.14782207 0.1112168 0.10489052
                0.3343035 0.3342455 0.2880437 0.15568811 0.1604730 0.06995111
                0.3066627 0.2807513 0.1494398 0.10650691 0.1370616 0.16258725
                0.2664826 0.2633140 0.1673342 0.12844107 0.1598561 0.16246321];
            
            D(2).Sigma_Y = 0.05*ones(16,6);
            
            D(2).condition = [1    0    0    0];
            
            D(3).t = [0  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30];
            
            D(3).Y = [...
                0.18118642 0.17117280 0.06921689 0.8339389 0.14309483 0.1119495
                0.15518812 0.14788162 0.38040487 0.9407151 0.18475791 0.7401733
                0.07325944 0.15237751 0.51660454 0.8499420 0.09395697 0.4317109
                0.13222423 0.14991635 0.69041908 0.9354031 0.11712634 0.4181825
                0.18434317 0.14962431 0.73298859 0.8944123 0.16154832 0.6943575
                0.13304805 0.12409090 0.84772786 0.9201726 0.21501215 0.2775180
                0.19558078 0.19270293 0.85919428 0.9246611 0.14410335 0.6985422
                0.16396667 0.18013749 0.81041944 0.8500798 0.15463602 0.4295047
                0.17767855 0.14559861 0.87204866 0.8830280 0.11152793 0.4746399
                0.12395204 0.12303433 0.92005853 0.9140457 0.15484416 0.5899543
                0.13706945 0.14579029 0.88166818 0.8935196 0.12884708 0.2104536
                0.13988815 0.09348882 0.84190569 0.8914391 0.07749091 0.7426955
                0.14811925 0.18067282 0.88561034 0.8828426 0.09152347 0.3317286
                0.13156674 0.18688392 0.86305399 0.9098012 0.06958192 0.5660655
                0.12220401 0.18730333 0.97629456 0.9202038 0.12613095 0.4910410
                0.14079056 0.12002588 0.88504812 0.8968963 0.20351941 0.3722563];
            
            D(3).Sigma_Y = 0.05*ones(16,6);
            
            D(3).condition = [0    1    0    0];
            
            D(4).t = [0  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30];
            
            D(4).Y = [...
                0.2045321 0.07560955 0.1640204 0.83166228 0.1457592 0.1879649
                0.6578856 0.72030634 0.3646157 0.62228680 0.1876390 0.7058527
                0.8727382 0.79903886 0.5831697 0.41246990 0.1924099 0.4796261
                0.8357646 0.88777519 0.7167020 0.29225199 0.3088524 0.3884216
                0.8144559 0.95385610 0.7399531 0.17231105 0.3646233 0.6760588
                0.7649645 0.85875447 0.7680542 0.18718742 0.4595235 0.2653407
                0.5937823 0.86985076 0.7847596 0.21682469 0.3860255 0.6389745
                0.6001083 0.71034331 0.8807288 0.17989607 0.4195041 0.4381725
                0.4909272 0.63099174 0.9054081 0.09502676 0.3755705 0.4303849
                0.4749780 0.56237922 0.9112906 0.15147265 0.3648891 0.6101374
                0.3513356 0.40745559 0.8654636 0.14254442 0.3041814 0.3246464
                0.2982424 0.32028891 0.8658242 0.09342384 0.3263700 0.6839370
                0.3013577 0.19068113 0.9092837 0.13758429 0.2392856 0.3283795
                0.3283447 0.20039130 0.9133164 0.11450109 0.2445276 0.5369872
                0.3221149 0.21700544 0.8690936 0.12761391 0.1700123 0.4998911
                0.3019580 0.28636982 0.8917507 0.09286719 0.1725707 0.3296038];
            
            D(4).Sigma_Y = 0.05*ones(16,6);
            
            D(4).condition = [1    1    0    0];
            
            D(5).t = [0  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30];
            
            D(5).Y = [...
                0.1632498 0.1479455 0.1575628 0.9343391 0.08689617 0.15362957
                0.6925302 0.8013656 0.2802410 0.8876928 0.11203316 0.19829660
                0.8492191 0.8804583 0.5279675 0.8020234 0.11350368 0.20716668
                0.8188087 0.8688254 0.6454666 0.8741729 0.15546884 0.13780278
                0.8075662 0.8962010 0.7192598 0.9433491 0.18198862 0.09098634
                0.6614125 0.8820678 0.8177839 0.8981930 0.15336392 0.13563796
                0.6122748 0.8159479 0.8128435 0.8540048 0.07683437 0.17942562
                0.5964614 0.7098129 0.8749133 0.8629443 0.19523545 0.10890770
                0.5203115 0.5989083 0.8561587 0.9122973 0.14203854 0.10378591
                0.4331483 0.5558090 0.7114718 0.8984247 0.11465383 0.13402169
                0.3735722 0.3486884 0.4803307 0.9206924 0.17063293 0.13564133
                0.3634435 0.3074528 0.3005396 0.9178086 0.16141023 0.16299146
                0.3163770 0.3145380 0.3053195 0.8925241 0.08437900 0.17649316
                0.2930702 0.2892375 0.2499740 0.9039381 0.15319574 0.12526669
                0.2859582 0.2580623 0.2360070 0.8561468 0.08830313 0.20297739
                0.2200915 0.2229491 0.1758068 0.8373987 0.14032307 0.08313549];
            
            D(5).Sigma_Y = 0.05*ones(16,6);
            
            D(5).condition = [1    0    1    0];
            
            D(6).t = [0  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30];
            
            D(6).Y = [...
                0.15802671 0.13643473 0.1800329 0.8484118 0.22663557 0.1555659
                0.14336998 0.13145709 0.3929189 0.8808722 0.19777232 0.7454059
                0.18040992 0.08567089 0.5400904 0.9302308 0.16982338 0.3367212
                0.12041775 0.14740326 0.6220083 0.8866694 0.13147906 0.3815005
                0.12697941 0.15359188 0.7508724 0.9446024 0.16607244 0.7251059
                0.13832002 0.20628729 0.7590274 0.9066425 0.17741022 0.2850733
                0.18236697 0.12795380 0.8237026 0.9295232 0.17331314 0.6564013
                0.10057642 0.07194705 0.9210295 0.9272505 0.14652365 0.3927551
                0.12600534 0.15459701 0.8368109 0.8935577 0.18190515 0.4273038
                0.17452037 0.09683811 0.9128952 0.9162230 0.15486223 0.5992842
                0.10454429 0.07742254 0.9598469 0.9137868 0.15418466 0.1964070
                0.13628964 0.21740637 0.8654660 0.8658182 0.20764707 0.6918639
                0.18508856 0.14711354 0.8812982 0.9109762 0.10824175 0.3930698
                0.11983420 0.20690178 0.8581597 0.9138032 0.03048901 0.5666084
                0.14001439 0.17982412 0.8702952 0.8788573 0.15310629 0.5137792
                0.05726728 0.11618360 0.8797293 0.8586403 0.13886133 0.2637185];
            
            D(6).Sigma_Y = 0.05*ones(16,6);
            
            D(6).condition = [0    1    1    0];
            
            D(7).t = [0  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30];
            
            D(7).Y = [...
                0.1354519 0.1535002 0.1355211 0.8948132 0.05830718 0.1331813
                0.6696247 0.8374504 0.3266970 0.9018727 0.04859938 0.6248577
                0.8569923 0.8792971 0.4591545 0.9034149 0.20142266 0.3991535
                0.8287373 0.8393217 0.6234086 0.8585909 0.35106319 0.3648269
                0.7883968 0.8529204 0.7264467 0.8813599 0.40183800 0.7142400
                0.7106160 0.8354582 0.7831088 0.9132600 0.37992788 0.2572520
                0.6793047 0.7827865 0.7864942 0.8869663 0.43525690 0.6229508
                0.5044554 0.7824801 0.8501095 0.9805354 0.42573238 0.3850916
                0.5465676 0.6065365 0.8549203 0.9267285 0.40746038 0.4725017
                0.4031440 0.5911519 0.8816015 0.8863326 0.35919635 0.5369462
                0.3644277 0.4043009 0.9506253 0.8743855 0.32692712 0.1892632
                0.3276202 0.2798091 0.8623968 0.8997609 0.27628588 0.7264940
                0.3531709 0.2515343 0.9306097 0.9151878 0.28073590 0.3068947
                0.3453740 0.2769636 0.9113310 0.8514286 0.23346398 0.5366783
                0.3526800 0.2103236 0.8655962 0.9037883 0.13719210 0.5167739
                0.3257124 0.2308489 0.9047394 0.8544476 0.15136109 0.2935793];
            
            D(7).Sigma_Y = 0.05*ones(16,6);
            
            D(7).condition = [1    1    1    0];
            
            D(8).t = [0  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30];
            
            D(8).Y = [...
                0.08280700 0.11553075 0.1472974 0.96285576 0.14497049 0.18884840
                0.17564872 0.12929936 0.3626649 0.62717539 0.13108944 0.19785061
                0.16561309 0.15084021 0.4754983 0.47963178 0.08875446 0.14219738
                0.11666629 0.16074090 0.6622509 0.34121863 0.11914046 0.19826344
                0.12026731 0.15630946 0.6436222 0.20653517 0.11793904 0.17812851
                0.15091108 0.18346866 0.8237203 0.22889624 0.09307591 0.20138569
                0.14032278 0.09836929 0.8282504 0.13764524 0.16710616 0.14555492
                0.11880250 0.16241808 0.8311246 0.24406630 0.14445627 0.14565816
                0.05929761 0.17015820 0.9341332 0.21564537 0.13859589 0.20944482
                0.17287625 0.21674687 0.8137794 0.08416098 0.13446667 0.16291694
                0.12925025 0.16079128 0.7916966 0.11918044 0.10265952 0.16558225
                0.21270642 0.12977357 0.8335809 0.09960676 0.16082558 0.17337739
                0.09927376 0.09580599 0.8863700 0.14684710 0.17718708 0.09018017
                0.10254795 0.15509262 0.8679626 0.11784491 0.13631146 0.16965269
                0.16297569 0.14581646 0.8885110 0.12995027 0.12279929 0.14879447
                0.06804257 0.17723955 0.9094174 0.14786865 0.14194530 0.17307476];
            
            D(8).Sigma_Y = 0.05*ones(16,6);
            
            D(8).condition = [1    0    0    1];
            
            D(9).t = [0  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30];
            
            D(9).Y = [...
                0.13316034 0.15323908 0.09803815 0.8527693 0.17100665 0.1257115
                0.11959537 0.16325565 0.34708573 0.9031208 0.08860702 0.6912611
                0.20986003 0.16529962 0.48657019 0.9013271 0.15110145 0.3889245
                0.10873692 0.14019480 0.65575738 0.8708600 0.09199577 0.3780692
                0.13048005 0.19814509 0.72354039 0.8910663 0.18477512 0.7900423
                0.16394083 0.11428609 0.76811811 0.8928696 0.16779222 0.3259523
                0.15878162 0.16321965 0.84381844 0.8870642 0.14254790 0.7397550
                0.10244079 0.08543513 0.84295148 0.9538362 0.17321911 0.4320170
                0.09371060 0.04413576 0.87647326 0.8901926 0.12877208 0.4130873
                0.12798497 0.13065782 0.88852442 0.9170711 0.08769015 0.6298504
                0.14273451 0.18143339 0.92679726 0.8408848 0.11502268 0.3149545
                0.13483210 0.14871571 0.91859347 0.9902228 0.15710659 0.8019711
                0.06203655 0.13298927 0.92963258 0.9329484 0.07242813 0.3026445
                0.15474661 0.13459618 0.92927695 0.9095623 0.13194513 0.5360921
                0.09798871 0.11136110 0.84310235 0.8947490 0.13654565 0.5563650
                0.09306310 0.14968652 0.86909037 0.9170778 0.13169712 0.3049163];
            
            D(9).Sigma_Y = 0.05*ones(16,6);
            
            D(9).condition = [0    1    0    1];
            
            D(10).t = [0  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30];
            
            D(10).Y = [...
                0.1143198 0.09963248 0.1555124 0.96880657 0.1055840 0.08355081
                0.1649738 0.17021499 0.2421683 0.55584174 0.1664282 0.70274948
                0.1545644 0.13505070 0.5307636 0.42850606 0.1760750 0.46763301
                0.1670297 0.06929900 0.6656633 0.25492390 0.3854150 0.39518933
                0.1719490 0.17362567 0.7707491 0.24598043 0.4247821 0.72558668
                0.1148706 0.16640169 0.8469331 0.22317895 0.5100184 0.28460000
                0.2386801 0.15844185 0.7899374 0.08670365 0.4730653 0.70365240
                0.1553100 0.13738471 0.8347681 0.19159311 0.5020870 0.41292336
                0.1882625 0.11955981 0.8743367 0.14719250 0.6386966 0.46428686
                0.1121880 0.09501891 0.8910863 0.14278949 0.7027587 0.57945952
                0.1456679 0.14302850 0.8442043 0.13812999 0.6978150 0.21135991
                0.1908594 0.14472561 0.8706587 0.07114821 0.7407568 0.72920072
                0.1033229 0.10535095 0.8787824 0.15172278 0.7779155 0.32307083
                0.1088544 0.18511048 0.8948015 0.12611032 0.7400662 0.56968804
                0.1689166 0.17694590 0.8542110 0.15519797 0.7646248 0.55091026
                0.1303701 0.11296231 0.8472156 0.19074312 0.8217805 0.34955475];
            
            D(10).Sigma_Y = 0.05*ones(16,6);
            
            D(10).condition = [1    1    0    1];
            
        case 5
            % Chen2009
            
            %A431
            
            % pAkt
            M = [... 
            0	0	0;...
            2.5	0	0;...
            5	24	0.26;...
            7.5	34	0.59;...
            10	40	2.2;...
            15	50	2;...
            30	64	1.4;...
            45	52	3.4;...
            60	49	2.2;...
            1.2e+02	57	7.3;...
            ];
        
        D(1).t = M(:,1)*60+1800; % minutes -> sec + preequi
        D(1).condition = [5e-9,0,0,0];
        D(1).Y(:,3) = M(:,2);
        D(1).Sigma_Y(:,3) = M(:,3);
        
        %pErbB1
        M= [...
            0	0	0;...
            2.5	65	7.3;...
            5	72	6.7;...
            7.5	75	8.8;...
            10	77	4.6;...
            15	83	4.6;...
            30	87	2.4;...
            45	98	2.3;...
            60	87	13;...
            1.2e+02	86	13;...
            ];
        
        D(1).Y(:,1) = M(:,2);
        D(1).Sigma_Y(:,1) = M(:,3);
        
        %pERK
        M = [...
            0	0	0;...
            2.5	60	0.36;...
            5	88	1.8;...
            7.5	97	1.4;...
            10	98	1.1;...
            15	86	3;...
            30	65	2.7;...
            45	56	2.8;...
            60	50	3.6;...
            1.2e+02	46	1.4;...
            ];
        
        D(1).Y(:,2) = M(:,2);
        D(1).Sigma_Y(:,2) = M(:,3);
        
        % pAkt
            M = [... 
                0	0	0;...
                2.5	45	4.6;...
                5	79	3.1;...
                7.5	86	4.5;...
                10	97	3.1;...
                15	91	2.8;...
                30	85	2.7;...
                45	72	2.5;...
                60	73	1.7;...
                1.2e+02	68	2.2;...
            ];
        
        D(2).t = M(:,1)*60+1800; % minutes -> sec + preequi
        D(2).condition = [0,0,0,5e-9];
        D(2).Y(:,3) = M(:,2);
        D(2).Sigma_Y(:,3) = M(:,3);
        
        %pErbB1
        M= [...
            0	0	0;...
            2.5	0	0;...
            5	0	0;...
            7.5	0	0;...
            10	0	0;...
            15	0	0;...
            30	0	0;...
            45	0	0;...
            60	0	0;...
            1.2e+02	0	0;...
            ];
        
        D(2).Y(:,1) = M(:,2);
        D(2).Sigma_Y(:,1) = M(:,3);
        
        %pERK
        M = [...
            0	0	0;...
            2.5	9.7	3.3;...
            5	46	6.7;...
            7.5	58	1.8;...
            10	54	3.9;...
            15	26	4.2;...
            30	27	2.6;...
            45	24	1.1;...
            60	23	0.84;...
            1.2e+02	17	2.2;...
            ];
        
        D(2).Y(:,2) = M(:,2);
        D(2).Sigma_Y(:,2) = M(:,3);
        
        % pAkt
        M = [...
            0	0	0;...
            2.5	1.7	0.88;...
            5	44	11;...
            7.5	43	5.5;...
            10	50	3;...
            15	45	4.4;...
            30	35	5;...
            45	40	14;...
            60	28	9.6;...
            1.2e+02	28	5.5;...
            ];
        
        D(3).t = M(:,1)*60+1800; % minutes -> sec + preequi
        D(3).condition = [1e-11,0,0,0];
        D(3).Y(:,3) = M(:,2);
        D(3).Sigma_Y(:,3) = M(:,3);
        
        %pErbB1
        M= [...
            0	0	0;...
            2.5	0	0;...
            5	0	0;...
            7.5	0	0;...
            10	0	0;...
            15	0	0;...
            30	0	0;...
            45	0	0;...
            60	0	0;...
            1.2e+02	0	0;...
            ];
        
        D(3).Y(:,1) = M(:,2);
        D(3).Sigma_Y(:,1) = M(:,3);
        
        %pERK
        M = [...
            0	0	0;...
            2.5	8.8	2.2;...
            5	43	8.8;...
            7.5	46	5.7;...
            10	40	3.2;...
            15	26	2.4;...
            30	27	0.28;...
            45	18	3.5;...
            60	11	2.4;...
            1.2e+02	2.7	0.65;...
            ];
        
        D(3).Y(:,2) = M(:,2);
        D(3).Sigma_Y(:,2) = M(:,3);
        
        % pAkt
        M = [...
            0	0	0;...
            2.5	1.7	0.88;...
            5	44	11;...
            7.5	43	5.5;...
            10	50	3;...
            15	45	4.4;...
            30	35	5;...
            45	40	14;...
            60	28	9.6;...
            1.2e+02	28	5.5;...
            ];
        
        D(4).t = M(:,1)*60+1800; % minutes -> sec + preequi
        D(4).condition = [0,0,0,1e-10];
        D(4).Y(:,3) = M(:,2);
        D(4).Sigma_Y(:,3) = M(:,3);
        
        %pErbB1
        M= [...
            0	0	0;...
            2.5	0	0;...
            5	0	0;...
            7.5	0	0;...
            10	0	0;...
            15	0	0;...
            30	0	0;...
            45	0	0;...
            60	0	0;...
            1.2e+02	0	0;...
            ];
        
        D(4).Y(:,1) = M(:,2);
        D(4).Sigma_Y(:,1) = M(:,3);
        
        %pERK
        M = [...
            0	0	0;...
            2.5	8.8	2.2;...
            5	43	8.8;...
            7.5	46	5.7;...
            10	40	3.2;...
            15	26	2.4;...
            30	27	0.28;...
            45	18	3.5;...
            60	11	2.4;...
            1.2e+02	2.7	0.65;...
            ];
        
        D(4).Y(:,2) = M(:,2);
        D(4).Sigma_Y(:,2) = M(:,3);
        
        for iD = 1:length(D)
            D(iD).Sigma_Y = D(iD).Sigma_Y + 0.1*ones(size(D(iD).Sigma_Y));
        end
       
        case 6
            % Kuhn2009
            load('./project/models/Kuhn2009_pnom.mat')
            load('./project/models/Kuhn2009_knom.mat')
            
            t = linspace(0,100,21);
            
            sol = simulate_Kuhn2009(t,log10(pnom),knom);
            
            D(1).t = t;
            D(1).condition = knom;
            D(1).Y = sol.y + 5*randn(size(sol.y));
            D(1).Sigma_Y = 0.1*ones(size(sol.y));
            
        case 7
            
            % Smith2013
            load('./project/models/Smith2013_pnom.mat')
            load('./project/models/Smith2013_knom.mat')
            
            t = linspace(0,10,21);
            
            sol = simulate_Smith2013(t,log10(pnom),knom);
            
            D(1).t = t;
            D(1).condition = knom;
            D(1).Y = sol.y + 1000*randn(size(sol.y));
            D(1).Sigma_Y = 1000*ones(size(sol.y));
        
    end
    
end


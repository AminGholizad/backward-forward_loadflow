branches = [%From	To	 R+jX (ohms)
    1	2	0.0922+1i*0.0470
    2	3	0.4930+1i*0.2511
    3	4	0.3660+1i*0.1864
    4	5	0.3811+1i*0.1941
    5	6	0.8190+1i*0.7070
    6	7	0.1872+1i*0.6188
    7	8	0.7114+1i*0.2351
    8	9	1.0300+1i*0.7400
    9	10	1.0440+1i*0.7400
    10	11	0.1966+1i*0.0650
    11	12	0.3744+1i*0.1298
    12	13	1.4680+1i*1.1550
    13	14	0.5416+1i*0.7129
    14	15	0.5910+1i*0.5260
    15	16	0.7463+1i*0.5450
    16	17	1.2890+1i*1.7210
    17	18	0.7320+1i*0.5740
    2	19	0.1640+1i*0.1565
    19	20	1.5042+1i*1.3554
    20	21	0.4095+1i*0.4784
    21	22	0.7089+1i*0.9373
    3	23	0.4512+1i*0.3083
    23	24	0.8980+1i*0.7091
    24	25	0.8960+1i*0.7011
    6	26	0.2030+1i*0.1034
    26	27	0.2842+1i*0.1447
    27	28	1.0590+1i*0.9337
    28	29	0.8042+1i*0.7006
    29	30	0.5075+1i*0.2585
    30	31	0.9744+1i*0.9630
    31	32	0.3105+1i*0.3619
    32	33	0.3410+1i*0.5302
    ];
spotloads = [%Bus Number	P+jQ(KVA)
    2	100+1i*60
    3	90+1i*40
    4	120+1i*80
    5	60+1i*30
    6	60+1i*20
    7	200+1i*100
    8	200+1i*100
    9	60+1i*20
    10	60+1i*20
    11	45+1i*30
    12	60+1i*35
    13	60+1i*35
    14	120+1i*80
    15	60+1i*10
    16	60+1i*20
    17	60+1i*20
    18	90+1i*40
    19	90+1i*40
    20	90+1i*40
    21	90+1i*40
    22	90+1i*40
    23	90+1i*50
    24	420+1i*200
    25	420+1i*200
    26	60+1i*25
    27	60+1i*25
    28	60+1i*20
    29	120+1i*70
    30	200+1i*600
    31	150+1i*70
    32	210+1i*100
    33	60+1i*40
    ];
spotloads=[spotloads(:,1) spotloads(:,2)*1e3];%P+jQ(VA)
Vb = 12.66;%(kV)
Vb = Vb*1e3;%(V)
rn = 1;% root node
[V, I] = LoadFlow(branches ,rn ,Vb , spotloads);
plot(abs(V(:,1)),abs(V(:,2)))
set(get(gca,'YAxis'),'Exponent',3)
ylabel('Voltage (V)')
xlabel('bus number')
xlim([1 33])
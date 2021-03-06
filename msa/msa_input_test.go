package msa

var testAlignedInput = []byte(`
>d1a1x__ b.63.1.1 (-)
PPDHLWVHQEGIYRDEYQRTWVAVVEE--E--T--SF---------LR----------ARVQQIQVPLG-------DAA
>gi|6678257|ref|NP_033363.1
HPNRLWIWEKHVYLDEFRRSWLPVVIK--S--N--EK---------FQ----------VILRQEDVTLG-------EAM
>gi|7305557|ref|NP_038800.1
PPRFLVCTRDDIYEDENGRQWVVAKVE--T--S--RSpygsrietcIT----------VHLQHMTTIPQ-------EPT
>gi|11415028|ref|NP_068801.1
HPDRLWAWEKFVYLDEKQHAWLPLTIEikD--R--LQ---------LR----------VLLRREDVVLG-------RPM
>gi|7305561|ref|NP_038804.1
----------GIYEDEHHRVWIAVNVE--T--S--HS---------SHgnrietcvt-VHLQHMTTLPQ-------EPT
>gi|7305553|ref|NP_038801.1
LPVYLVSVRLGIYEDEHHRVWIVANVE--TshS--SH---------GN----------RRRTHVTVHLW-------KLI
>gi|27668591|ref|XP_234504.1
-PDRLWLWEKHVYLDEFRRSWLPIVIK--S--N--GK---------FQ----------VIMRQKDVILG-------DSM
`)

var testBadAlignedInput = []byte(`
>d1a1x__ b.63.1.1 (-)
PPLW
>gi|6678257|ref|NP_033363.1
HPNRLWIWEKHVYLDEFRRSWLPVVIK--S--N--EK---------FQ----------VILRQEDVTLG-------EAM
>gi|7305557|ref|NP_038800.1
PPRFLVCTRDDIYEDENGRQWVVAKVE--T--S--RSpygsrietcIT----------VHLQHMTTIPQ-------EPT
`)

var inpAlignFasta = []string{
	">1",
	"PPDHLWVHQEGIYRDEYQRTWVAVVEE--E--T--SF---------LR----------ARVQQIQVPLG---" +
		"----DAARPSHLLTS-----QL",
	">2",
	"HPNRLWIWEKHVYLDEFRRSWLPVVIK--S--N--EK---------FQ----------VILRQEDVTLG---" +
		"----EAMSPSQLVPY-----EL",
	">3",
	"PPRFLVCTRDDIYEDENGRQWVVAKVE--T--S--RSpygsrietcIT----------VHLQHMTTIPQ---" +
		"----EPTPQQPINNN-----SL",
	">4",
	"HPDRLWAWEKFVYLDEKQHAWLPLTIEikD--R--LQ---------LR----------VLLRREDVVLG---" +
		"----RPMTPTQIGPS-----LL",
	">5",
	"----------GIYEDEHHRVWIAVNVE--T--S--HS---------SHgnrietcvt-VHLQHMTTLPQ---" +
		"----EPTPQQPINNN-----SL",
	">6",
	"LPVYLVSVRLGIYEDEHHRVWIVANVE--TshS--SH---------GN----------RRRTHVTVHLW---" +
		"----KLIPQQVIPFNplnydFL",
	">7",
	"-PDRLWLWEKHVYLDEFRRSWLPIVIK--S--N--GK---------FQ----------VIMRQKDVILG---" +
		"----DSMTPSQLVPY-----EL",
	">8",
	"-PHILTLRTHGIYEDEHHRLWVVLDLQ--A--ShlSF---------SN----------RLLIYLTVYLQqgv" +
		"afplESTPPSPMNLN-----GL",
	">9",
	"PPCFLVCTRDDIYEDEHGRQWVAAKVE--T--S--SH---------SPycskietcvtVHLWQMTTLFQ---" +
		"----EPSPDSLKTFN-----FL",
	">10",
	"---------PGFYEDEHHRLWMVAKLE--T--C--SH---------SPycnkietcvtVHLWQMTRYPQ---" +
		"----EPAPYNPMNYN-----FL",
	">11",
	"---A--------------------------------------------------------------------" +
		"---------------------L",
}

var inpAlignA2M = []string{
	">1",
	"PPDHLWVHQEGIYRDEYQRTWVAVVEE..E..T..SF.........LR..........ARVQQIQVPLG..." +
		"....DAARPSHLLTS.....QL",
	">2",
	"HPNRLWIWEKHVYLDEFRRSWLPVVIK..S..N..EK.........FQ..........VILRQEDVTLG..." +
		"....EAMSPSQLVPY.....EL",
	">3",
	"PPRFLVCTRDDIYEDENGRQWVVAKVE..T..S..RSpygsrietcIT..........VHLQHMTTIPQ..." +
		"....EPTPQQPINNN.....SL",
	">4",
	"HPDRLWAWEKFVYLDEKQHAWLPLTIEikD..R..LQ.........LR..........VLLRREDVVLG..." +
		"....RPMTPTQIGPS.....LL",
	">5",
	"----------GIYEDEHHRVWIAVNVE..T..S..HS.........SHgnrietcvt.VHLQHMTTLPQ..." +
		"....EPTPQQPINNN.....SL",
	">6",
	"LPVYLVSVRLGIYEDEHHRVWIVANVE..TshS..SH.........GN..........RRRTHVTVHLW..." +
		"....KLIPQQVIPFNplnydFL",
	">7",
	"-PDRLWLWEKHVYLDEFRRSWLPIVIK..S..N..GK.........FQ..........VIMRQKDVILG..." +
		"....DSMTPSQLVPY.....EL",
	">8",
	"-PHILTLRTHGIYEDEHHRLWVVLDLQ..A..ShlSF.........SN..........RLLIYLTVYLQqgv" +
		"afplESTPPSPMNLN.....GL",
	">9",
	"PPCFLVCTRDDIYEDEHGRQWVAAKVE..T..S..SH.........SPycskietcvtVHLWQMTTLFQ..." +
		"....EPSPDSLKTFN.....FL",
	">10",
	"---------PGFYEDEHHRLWMVAKLE..T..C..SH.........SPycnkietcvtVHLWQMTRYPQ..." +
		"....EPAPYNPMNYN.....FL",
	">11",
	"---A-----------------------..-..-..--.........--..........-----------..." +
		"....-----------.....-L",
}

var inpAlignA3M = []string{
	">1",
	"PPDHLWVHQEGIYRDEYQRTWVAVVEEETSFLRARVQQIQVPLGDAARPSHLLTSQL",
	">2",
	"HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL",
	">3",
	"PPRFLVCTRDDIYEDENGRQWVVAKVETSRSpygsrietcITVHLQHMTTIPQEPTPQQPINNNSL",
	">4",
	"HPDRLWAWEKFVYLDEKQHAWLPLTIEikDRLQLRVLLRREDVVLGRPMTPTQIGPSLL",
	">5",
	"----------GIYEDEHHRVWIAVNVETSHSSHgnrietcvtVHLQHMTTLPQEPTPQQPINNNSL",
	">6",
	"LPVYLVSVRLGIYEDEHHRVWIVANVETshSSHGNRRRTHVTVHLWKLIPQQVIPFNplnydFL",
	">7",
	"-PDRLWLWEKHVYLDEFRRSWLPIVIKSNGKFQVIMRQKDVILGDSMTPSQLVPYEL",
	">8",
	"-PHILTLRTHGIYEDEHHRLWVVLDLQAShlSFSNRLLIYLTVYLQqgvafplESTPPSPMNLNGL",
	">9",
	"PPCFLVCTRDDIYEDEHGRQWVAAKVETSSHSPycskietcvtVHLWQMTTLFQEPSPDSLKTFNFL",
	">10",
	"---------PGFYEDEHHRLWMVAKLETCSHSPycnkietcvtVHLWQMTRYPQEPAPYNPMNYNFL",
	">11",
	"---A----------------------------------------------------L",
}

var alignFasta = []string{
	"PPDHLWVHQEGIYRDEYQRTWVAVVEE--E--T--SF---------LR----------ARVQQIQVPLG---" +
		"----DAARPSHLLTS-----QL",
	"HPNRLWIWEKHVYLDEFRRSWLPVVIK--S--N--EK---------FQ----------VILRQEDVTLG---" +
		"----EAMSPSQLVPY-----EL",
	"PPRFLVCTRDDIYEDENGRQWVVAKVE--T--S--RSpygsrietcIT----------VHLQHMTTIPQ---" +
		"----EPTPQQPINNN-----SL",
	"HPDRLWAWEKFVYLDEKQHAWLPLTIEikD--R--LQ---------LR----------VLLRREDVVLG---" +
		"----RPMTPTQIGPS-----LL",
	"----------GIYEDEHHRVWIAVNVE--T--S--HS---------SHgnrietcvt-VHLQHMTTLPQ---" +
		"----EPTPQQPINNN-----SL",
	"LPVYLVSVRLGIYEDEHHRVWIVANVE--TshS--SH---------GN----------RRRTHVTVHLW---" +
		"----KLIPQQVIPFNplnydFL",
	"-PDRLWLWEKHVYLDEFRRSWLPIVIK--S--N--GK---------FQ----------VIMRQKDVILG---" +
		"----DSMTPSQLVPY-----EL",
	"-PHILTLRTHGIYEDEHHRLWVVLDLQ--A--ShlSF---------SN----------RLLIYLTVYLQqgv" +
		"afplESTPPSPMNLN-----GL",
	"PPCFLVCTRDDIYEDEHGRQWVAAKVE--T--S--SH---------SPycskietcvtVHLWQMTTLFQ---" +
		"----EPSPDSLKTFN-----FL",
	"---------PGFYEDEHHRLWMVAKLE--T--C--SH---------SPycnkietcvtVHLWQMTRYPQ---" +
		"----EPAPYNPMNYN-----FL",
	"---A--------------------------------------------------------------------" +
		"---------------------L",
}

var alignA2M = []string{
	"PPDHLWVHQEGIYRDEYQRTWVAVVEE..E..T..SF.........LR..........ARVQQIQVPLG..." +
		"....DAARPSHLLTS.....QL",
	"HPNRLWIWEKHVYLDEFRRSWLPVVIK..S..N..EK.........FQ..........VILRQEDVTLG..." +
		"....EAMSPSQLVPY.....EL",
	"PPRFLVCTRDDIYEDENGRQWVVAKVE..T..S..RSpygsrietcIT..........VHLQHMTTIPQ..." +
		"....EPTPQQPINNN.....SL",
	"HPDRLWAWEKFVYLDEKQHAWLPLTIEikD..R..LQ.........LR..........VLLRREDVVLG..." +
		"....RPMTPTQIGPS.....LL",
	"----------GIYEDEHHRVWIAVNVE..T..S..HS.........SHgnrietcvt.VHLQHMTTLPQ..." +
		"....EPTPQQPINNN.....SL",
	"LPVYLVSVRLGIYEDEHHRVWIVANVE..TshS..SH.........GN..........RRRTHVTVHLW..." +
		"....KLIPQQVIPFNplnydFL",
	"-PDRLWLWEKHVYLDEFRRSWLPIVIK..S..N..GK.........FQ..........VIMRQKDVILG..." +
		"....DSMTPSQLVPY.....EL",
	"-PHILTLRTHGIYEDEHHRLWVVLDLQ..A..ShlSF.........SN..........RLLIYLTVYLQqgv" +
		"afplESTPPSPMNLN.....GL",
	"PPCFLVCTRDDIYEDEHGRQWVAAKVE..T..S..SH.........SPycskietcvtVHLWQMTTLFQ..." +
		"....EPSPDSLKTFN.....FL",
	"---------PGFYEDEHHRLWMVAKLE..T..C..SH.........SPycnkietcvtVHLWQMTRYPQ..." +
		"....EPAPYNPMNYN.....FL",
	"---A-----------------------..-..-..--.........--..........-----------..." +
		"....-----------.....-L",
}

var alignA3M = []string{
	"PPDHLWVHQEGIYRDEYQRTWVAVVEEETSFLRARVQQIQVPLGDAARPSHLLTSQL",
	"HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL",
	"PPRFLVCTRDDIYEDENGRQWVVAKVETSRSpygsrietcITVHLQHMTTIPQEPTPQQPINNNSL",
	"HPDRLWAWEKFVYLDEKQHAWLPLTIEikDRLQLRVLLRREDVVLGRPMTPTQIGPSLL",
	"----------GIYEDEHHRVWIAVNVETSHSSHgnrietcvtVHLQHMTTLPQEPTPQQPINNNSL",
	"LPVYLVSVRLGIYEDEHHRVWIVANVETshSSHGNRRRTHVTVHLWKLIPQQVIPFNplnydFL",
	"-PDRLWLWEKHVYLDEFRRSWLPIVIKSNGKFQVIMRQKDVILGDSMTPSQLVPYEL",
	"-PHILTLRTHGIYEDEHHRLWVVLDLQAShlSFSNRLLIYLTVYLQqgvafplESTPPSPMNLNGL",
	"PPCFLVCTRDDIYEDEHGRQWVAAKVETSSHSPycskietcvtVHLWQMTTLFQEPSPDSLKTFNFL",
	"---------PGFYEDEHHRLWMVAKLETCSHSPycnkietcvtVHLWQMTRYPQEPAPYNPMNYNFL",
	"---A----------------------------------------------------L",
}

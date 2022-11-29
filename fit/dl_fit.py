import ROOT
import time

def jpsi_dl_fit():
        
    file = ROOT.TFile.Open("data_root_files/RunB_HLT_Dimuon25.root")

    ## Decay lenght fit

    # Decay lenght 
    jpsi_dl = ROOT.RooRealVar("jpsi_dl", "", -0.5, 0.8)
    #jpsi_dl = ROOT.RooRealVar("jpsi_dl", "", 0.05, 0.6)

    jpsi_dlErr = ROOT.RooRealVar("jpsi_dlErr", "", 4 ,0.001, 10)
    #jpsi_dlErr.setConstant(True)

    root_data = file.asso

    # Data
    data = ROOT.RooDataSet("data", "", ROOT.RooArgSet(jpsi_dl,), ROOT.RooFit.Import(root_data))

    ## First try: Gaussian + exponential
    #sigma_dl = ROOT.RooRealVar("sigma_dl", "", 0.01, 5)

    ## PDFs for prompt signal: Resolution function 

    # Prompt Resolution 
    mean_prompt = ROOT.RooRealVar("mean_prompt", "", 0.0005, -0.1, 0.1)
    sigma_prompt = ROOT.RooRealVar("sigma_prompt", "", 0.0060, 0, 0.1)
    frac_prompt = ROOT.RooRealVar("frac_prompt", "", 0.3, 0, 1.0)
    resolution_prompt = ROOT.RooGaussModel("resolution_prompt", "", jpsi_dl, mean_prompt, sigma_prompt, jpsi_dlErr)

    # Non Prompt Resolution 
    mean_non_prompt = ROOT.RooRealVar("mean_non_prompt", "", 0.0005, -1, 1)
    sigma_non_prompt = ROOT.RooRealVar("sigma_non_prompt", "", 0.0040, 0, 0.1)
    frac_non_prompt = ROOT.RooRealVar("frac_non_prompt", "", 0.88, 0, 1.0)
    resolution_non_prompt = ROOT.RooGaussModel("resolution_non_prompt", "", jpsi_dl, mean_non_prompt, sigma_non_prompt, jpsi_dlErr)

    # Gaussian 
    mean_pv = ROOT.RooRealVar("mean_pv", "", 0.00005, -0.3, 0.3)
    sigma_pv = ROOT.RooRealVar("sigma_pv", "", 0.0030, 0, 0.3)
    gauss_pv = ROOT.RooGaussian("gauss_PV", "", jpsi_dl, mean_pv, sigma_pv)

    # Prompt signal: Resolution function
    prompt_signal = ROOT.RooAddPdf("prompt_signal", "", ROOT.RooArgList(resolution_prompt, gauss_pv),
                        ROOT.RooArgList(frac_prompt), ROOT.kTRUE)

    ## PDFs for nomprompt signal

    ## Reference to add model: https://root.cern/doc/master/rf209__anaconv_8C.html

    # Exponential decay
    #tau = ROOT.RooRealVar("tau", "", 0.05, -10, 2)
    tau = ROOT.RooRealVar("tau", "", 0.065, -10, 2)
    #tau.setConstant(True)

    # PDF for nompromt signal
    frac_exp_1 = ROOT.RooRealVar("frac_exp_1", "", 0.0001, 0, 1)
    exp_decay = ROOT.RooDecay("exp_decay_sing", "", jpsi_dl, tau, resolution_non_prompt, ROOT.RooDecay.SingleSided)

    #c_exp = ROOT.RooRealVar("exp_decay_sing","Exp arg constant", -10, 0) 
    #exp_decay = ROOT.RooExponential("exp_decay_sing","BkgPDF", jpsi_dl, c_exp)
    #exp_decay = ROOT.RooExponential("exp_decay", "", jpsi_dl, tau)
    """ exp_decay_doub = ROOT.RooDecay("exp_decay_doub", "", jpsi_dl, tau, resolution_non_prompt, ROOT.RooDecay.DoubleSided)
    exp_decay_flip = ROOT.RooDecay("exp_decay_flip", "", jpsi_dl, tau, resolution_non_prompt, ROOT.RooDecay.Flipped)
    exp_decay = ROOT.RooAddPdf("exp_decay", "", ROOT.RooArgList(exp_decay_sing, exp_decay_doub, exp_decay_flip),
                        ROOT.RooArgList(frac_exp_1, frac_exp_2), ROOT.kTRUE) """

    # Model pdf
    model = ROOT.RooAddPdf("model", "", ROOT.RooArgList(prompt_signal, exp_decay),
                                        ROOT.RooArgList(frac_non_prompt), ROOT.kTRUE)

    result = model.fitTo(data, ROOT.RooFit.Save(), ROOT.RooFit.ConditionalObservables(jpsi_dlErr))

    frame1 = jpsi_dl.frame(ROOT.RooFit.Title(""))
    frame1.GetXaxis().SetTitle("#l_{J/\psi}\ [mm]")

    data.plotOn(frame1, ROOT.RooFit.Name("jpsi_dl"), ROOT.RooFit.Name("Data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

    # Colors and styles
    colors = {"model" : 2, "prompt" : 4, "non-prompt" : 3, "background": 1}
    styles = {"model" : 1, "prompt" : 1, "non-prompt" : 2, "background": 2}

    # Prompt
    prompt_signal.plotOn(frame1, ROOT.RooFit.Name("prompt_signal"), ROOT.RooFit.LineStyle(styles["prompt"]),
                    ROOT.RooFit.LineColor(colors["prompt"]))

    # Non-Prompt
    exp_decay.plotOn(frame1, ROOT.RooFit.Name("non_prompt_signal"), ROOT.RooFit.LineStyle(styles["non-prompt"]),
                    ROOT.RooFit.LineColor(colors["non-prompt"]))

    # Model
    model.plotOn(frame1, ROOT.RooFit.Name("model"), ROOT.RooFit.ProjWData(data), ROOT.RooFit.LineStyle(styles["model"]),
                    ROOT.RooFit.LineColor(colors["model"])) 

    can = ROOT.TCanvas("can", "histograms", 850, 600)
    can.Draw()
    leg_dl = ROOT.TLegend(0.7, 0.7, 0.88, 0.89)
    leg_dl.AddEntry(frame1.findObject("Data"), "Data", "LEP")
    leg_dl.AddEntry(frame1.findObject("model"), "Model Fit", "L")
    leg_dl.AddEntry(frame1.findObject("prompt_signal"), "Prompt Fit", "L")
    leg_dl.AddEntry(frame1.findObject("non_prompt_signal"), "Non-Prompt fit", "L")

    frame1.Draw()
    leg_dl.Draw("same")

    can.SetLogy()
    can.SaveAs("jpsi_dl_alone.png")

    ## Jpsi dl pull distribution
    cpjdl = ROOT.TCanvas("Jpsi decay length pull canvas")

    # Creates the pull histogram
    histpull_jpsi_dl = frame1.pullHist()

    # New frame to draw pull distribution
    frame_pull_jpsi_dl = jpsi_dl.frame(ROOT.RooFit.Title("Pull Distribution"))
    frame_pull_jpsi_dl.GetXaxis().SetTitle("#l_{J/\psi}\ [mm]")

    # Add the distribution to the frame
    frame_pull_jpsi_dl.addPlotable(histpull_jpsi_dl, "P")

    ROOT.gPad.SetLeftMargin(0.15)
    frame_pull_jpsi_dl.GetYaxis().SetTitleOffset(1.6)
    frame_pull_jpsi_dl.Draw()

    right = ROOT.TLatex()
    right.SetNDC()
    right.SetTextFont(43)
    right.SetTextSize(30)
    right.SetTextAlign(13)
    right.DrawLatex(.10,.95,"#bf{CMS} #scale[0.7]{#it{Preliminary}}")
    right.SetTextSize(20)
    right.DrawLatex(.80,.94 , "4.80 fb^{-1}")

    cpjdl.SaveAs('pull_jpsi_dl_alone.png')
    
    # Jpsi dl chi square  
    nparm_jpsi_dl = 8
    chi_square_jpsi_dl = frame1.chiSquare(nparm_jpsi_dl)
    print(f"Xi square for J/psi decay length is: {chi_square_jpsi_dl}")


jpsi_dl_fit()















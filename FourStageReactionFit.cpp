//==============================================================================================================================
//===================================          4-Stage multiple fit analysis     ===============================================
// This script do a simultaneous fit analysis of a four stage kinetic scheme. It has been developed to extract the kinetic
// parameters reported in the article entitled: Determination of kinetic properties in unimolecular dissociation of
// complex systems from graph-theory based analysis of an ensemble of reactive trajectories. Please for further details 
// please read the  article.
//
//  The installation of ROOT can be done easily through the following steps: https://iscinumpy.gitlab.io/post/root-conda/
//
// Author: Ariel Francis Perez Mellor
// Last Update: 2021-05-22
//
// input files:   00prob-SP        =  Probability of the starting point state as a function of time. <time> <Prob> <ePprob> <#events>     
//                01prob-INT       =  Probability of the intermediate state as a function of time.
//                03prob-M03-total =  Probability of the state A as a function of time.
//                04prob-PF-M03    =  Probability of the state B as a function of time.
//
//   parameters
//   p[0] = kis 
//   p[1] = ksi        
//   p[2] = kas     
//   p[3] = kai    
//   p[4] = kas + kbs 
//   p[5] = kai + kbi 
//   p[6] = t0   (it is not included in the article, it is a time shift ) 
//   p[7] = kba
//
// ouput files: MULTIFIT-RESULTS_FIT    = display all the rate constants values as well as the statistical errors 
//              MULTIFIT-RESULTS_LINES  = display the population of the states obtained from the fit
//              MultipleFit.pdf         = display the figures resulting from the fit    
//
//==============================================================================================================================
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include <TStopwatch.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TVectorT.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TPave.h"

#include "FourFunctions.h"


using namespace std;

//  here we define the shared parameters
    int iparCHAN[8] = {                        
        0,        //    kis 
        1,        //    ksi 
        2,        //    kas  
        3,        //    kai 
        4,        //    kas + kbs 
        5,        //    kai + kbi  
        6,        //    t0   
        7,        //    kba
        };    
    int iparREST[8] = {                        
        0,        //    kis 
        1,        //    ksi 
        2,        //    kas  
        3,        //    kai 
        4,        //    kas + kbs 
        5,        //    kai + kbi  
        6,        //    t0   
        7,        //    kba
        };    
    int iparINT[8] = {                        
        0,        //    kis 
        1,        //    ksi 
        2,        //    kas  
        3,        //    kai 
        4,        //    kas + kbs 
        5,        //    kai + kbi  
        6,        //    t0   
        7,        //    kba
        };    
    int iparNO_REACT[8] = {                        
        0,        //    kis 
        1,        //    ksi 
        2,        //    kas  
        3,        //    kai 
        4,        //    kas + kbs 
        5,        //    kai + kbi  
        6,        //    t0   
        7,        //    kba
        };                

//  here we construct the global chi structure
    struct GlobalChi2 { GlobalChi2(  ROOT::Math::IMultiGenFunction & f0, ROOT::Math::IMultiGenFunction & f1, ROOT::Math::IMultiGenFunction & f2, ROOT::Math::IMultiGenFunction & f3) 
                                    : fChi2_0(&f0), fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}
       double operator() (const double *par) const {
            
            double p0[8];
            for (int i = 0; i < 8; ++i) p0[i] = par[ iparCHAN[i] ];
        
            double p1[8];
            for (int i = 0; i < 8; ++i) p1[i] = par[ iparREST[i] ];
            
            double p2[8];
            for (int i = 0; i < 8; ++i) p2[i] = par[ iparINT[i] ];
            
            double p3[8];
            for (int i = 0; i < 8; ++i) p3[i] = par[ iparNO_REACT[i] ];
                    
            return (*fChi2_0)(p0) + (*fChi2_1)(p1) + (*fChi2_2)(p2)+(*fChi2_3)(p3) ;
        }

        const  ROOT::Math::IMultiGenFunction * fChi2_0;
        const  ROOT::Math::IMultiGenFunction * fChi2_1;
        const  ROOT::Math::IMultiGenFunction * fChi2_2;
        const  ROOT::Math::IMultiGenFunction * fChi2_3;                  
        };


void FourStageReactionFit()
{

    TStopwatch timer;
    timer.Start();    

// loading input data 
    char CHAN[]="03prob-M03-total";
    double cTIME, PROB, ePROB, tmp;
    TVectorD PROB_CHAN, ePROB_CHAN, TIME, eTIME;
    ifstream fin_CHAN(CHAN);

    int j=0;
    while(fin_CHAN >> cTIME >> PROB >> ePROB >> tmp)
    {
            PROB_CHAN.ResizeTo(j+1);
            ePROB_CHAN.ResizeTo(j+1);
            TIME.ResizeTo(j+1);
            eTIME.ResizeTo(j+1);
                        
            PROB_CHAN[j]=PROB;
            ePROB_CHAN[j]=ePROB;
            TIME[j]=cTIME;
            eTIME[j]=0.0;
            
            j++;
    }
    fin_CHAN.close();    fin_CHAN.clear();    
    
    char REST[]="04prob-PF-M03";
    TVectorD PROB_REST, ePROB_REST;
    ifstream fin_REST(REST);

    j=0;
    while(fin_REST >> cTIME >> PROB >> ePROB >> tmp)
    {
            PROB_REST.ResizeTo(j+1);
            ePROB_REST.ResizeTo(j+1);
                        
            PROB_REST[j]=PROB;
            ePROB_REST[j]=ePROB;
            
            j++;
    }
    fin_REST.close();    fin_REST.clear();    

    char INT[]="01prob-INT";
    TVectorD PROB_INT, ePROB_INT;
    ifstream fin_INT(INT);

    j=0;
    while(fin_INT  >> cTIME >> PROB >> ePROB >> tmp)
    {
            PROB_INT.ResizeTo(j+1);
            ePROB_INT.ResizeTo(j+1);
                        
            PROB_INT[j]=PROB;
            ePROB_INT[j]=ePROB;
        
            j++;
    }
    fin_INT.close();    fin_INT.clear();    

    char NO_REACT[]="00prob-SP";
    TVectorD PROB_NO_REACT, ePROB_NO_REACT;
    ifstream fin_NO_REACT(NO_REACT);    
    
    j=0;
    while(fin_NO_REACT  >> cTIME >> PROB >> ePROB >> tmp)
    {
            PROB_NO_REACT.ResizeTo(j+1);
            ePROB_NO_REACT.ResizeTo(j+1);
                        
            PROB_NO_REACT[j]=PROB;
            ePROB_NO_REACT[j]=ePROB;
        
            j++;
    }
    fin_NO_REACT.close();    fin_NO_REACT.clear();    

//  starting the fitting
    gROOT->SetStyle("Plain");
                
        if(gROOT->FindObject("DATA_PROB_CHAN")!=0)
        {
           gROOT->ProcessLine("delete DATA_PROB_CHAN");
        }
        
        TGraphErrors *DATA_PROB_CHAN=new TGraphErrors( TIME, PROB_CHAN, eTIME, ePROB_CHAN);

        DATA_PROB_CHAN->SetLineColor(1);
        DATA_PROB_CHAN->SetMarkerStyle(20);
        DATA_PROB_CHAN->SetMarkerColor(1);
        DATA_PROB_CHAN->SetMarkerSize(0.8);
//______________________________________________________________________________________________    
     gROOT->SetStyle("Plain");
             
        if(gROOT->FindObject("DATA_PROB_REST")!=0)
        {
           gROOT->ProcessLine("delete DATA_PROB_REST");
        }
        
        TGraphErrors *DATA_PROB_REST=new TGraphErrors( TIME, PROB_REST, eTIME, ePROB_REST);

        DATA_PROB_REST->SetLineColor(1);
        DATA_PROB_REST->SetMarkerStyle(20);
        DATA_PROB_REST->SetMarkerColor(1);
        DATA_PROB_REST->SetMarkerSize(0.8);
//______________________________________________________________________________________________    
        if(gROOT->FindObject("DATA_PROB_INT")!=0)
        {
           gROOT->ProcessLine("delete DATA_PROB_INT");
        }
        
        TGraphErrors *DATA_PROB_INT=new TGraphErrors( TIME, PROB_INT, eTIME, ePROB_INT);

        DATA_PROB_INT->SetLineColor(1);
        DATA_PROB_INT->SetMarkerStyle(20);
        DATA_PROB_INT->SetMarkerColor(1);
        DATA_PROB_INT->SetMarkerSize(0.8);
//______________________________________________________________________________________________  
        if(gROOT->FindObject("DATA_PROB_NO_REACT")!=0)
        {
           gROOT->ProcessLine("delete DATA_PROB_NO_REACT");
        }        
        
        TGraphErrors *DATA_PROB_NO_REACT=new TGraphErrors( TIME, PROB_NO_REACT, eTIME, ePROB_NO_REACT);

        DATA_PROB_NO_REACT->SetLineColor(1);
        DATA_PROB_NO_REACT->SetMarkerStyle(20);
        DATA_PROB_NO_REACT->SetMarkerColor(1);
        DATA_PROB_NO_REACT->SetMarkerSize(0.8);
//______________________________________________________________________________________________   
        if(gROOT->FindObject("FIT_CHAN")!=0)
        {
           gROOT->ProcessLine("delete FIT_CHAN");
        }
        TF1 *FIT_CHAN=new TF1("FIT_CHAN", ADJ_CHAN_4th, TIME.Min(), TIME.Max(), 8);
//______________________________________________________________________________________________    
        if(gROOT->FindObject("FIT_REST")!=0)
        {
           gROOT->ProcessLine("delete FIT_REST");
        }
        TF1 *FIT_REST=new TF1("FIT_REST", ADJ_REST_4th, TIME.Min(), TIME.Max(), 8);
//______________________________________________________________________________________________    
        if(gROOT->FindObject("FIT_INT")!=0)
        {
           gROOT->ProcessLine("delete FIT_INT");
        }
        
        TF1 *FIT_INT=new TF1("FIT_INT", ADJ_INT_4th, TIME.Min(), TIME.Max(), 8);
//______________________________________________________________________________________________
        if(gROOT->FindObject("FIT_NO_REACT")!=0)
        {
           gROOT->ProcessLine("delete FIT_NO_REACT");
        }
        
        TF1 *FIT_NO_REACT=new TF1("FIT_NO_REACT", ADJ_NO_REACT_4th, TIME.Min(), TIME.Max(), 8);    
//______________________________________________________________________________________________     
        ROOT::Math::WrappedMultiTF1 wFIT_CHAN(*FIT_CHAN, 1);
        ROOT::Math::WrappedMultiTF1 wFIT_REST(*FIT_REST, 1);
        ROOT::Math::WrappedMultiTF1 wFIT_INT(*FIT_INT, 1);
        ROOT::Math::WrappedMultiTF1 wFIT_NO_REACT(*FIT_NO_REACT, 1);
    
        ROOT::Fit::DataOptions opt;

        ROOT::Fit::DataRange rangeFIT_CHAN;
        rangeFIT_CHAN.SetRange(TIME.Min(), TIME.Max());
        ROOT::Fit::BinData dataFIT_CHAN(opt, rangeFIT_CHAN);
        ROOT::Fit::FillData(dataFIT_CHAN, DATA_PROB_CHAN);
//______________________________________________________________________________________________    
        ROOT::Fit::DataRange rangeFIT_REST;
        rangeFIT_REST.SetRange(TIME.Min(), TIME.Max());
        ROOT::Fit::BinData dataFIT_REST(opt, rangeFIT_REST);
        ROOT::Fit::FillData(dataFIT_REST, DATA_PROB_REST);
//______________________________________________________________________________________________    
        ROOT::Fit::DataRange rangeFIT_INT;
        rangeFIT_INT.SetRange(TIME.Min(), TIME.Max());
        ROOT::Fit::BinData dataFIT_INT(opt, rangeFIT_INT);
        ROOT::Fit::FillData(dataFIT_INT, DATA_PROB_INT);
//______________________________________________________________________________________________    
        ROOT::Fit::DataRange rangeFIT_NO_REACT;
        rangeFIT_NO_REACT.SetRange(TIME.Min(), TIME.Max());
        ROOT::Fit::BinData dataFIT_NO_REACT(opt, rangeFIT_NO_REACT);
        ROOT::Fit::FillData(dataFIT_NO_REACT, DATA_PROB_NO_REACT);
//______________________________________________________________________________________________        
        ROOT::Fit::Chi2Function chi2_CHAN(dataFIT_CHAN, wFIT_CHAN);
        ROOT::Fit::Chi2Function chi2_REST(dataFIT_REST, wFIT_REST);
        ROOT::Fit::Chi2Function chi2_INT(dataFIT_INT, wFIT_INT);
        ROOT::Fit::Chi2Function chi2_NO_REACT(dataFIT_NO_REACT, wFIT_NO_REACT);
        
        GlobalChi2 globalChi2(chi2_CHAN, chi2_REST, chi2_INT, chi2_NO_REACT);

        ROOT::Fit::Fitter fitter;
        
        const int Npar =8;
        double par0[Npar] = {    
        
        //   initializing the parameters            
        0.000246958,              // kis         0
        2.04081e-05,              // ksi         1
        0,                        // kas         2
        1.3888e-05,               // kai         3
        2.08424e-05,              // kas + kbs   4
        0.000169047,              // kai + kbi   5
        0,                        // t0          6
        0,                        // kba         7
        
        };        

//  parameter settings 
    fitter.Config().SetParamsSettings(8, par0);
  
    fitter.Config().ParSettings(0).Fix();
    fitter.Config().ParSettings(1).Fix();
    fitter.Config().ParSettings(2).Fix();
//  fitter.Config().ParSettings(3).Fix();
    fitter.Config().ParSettings(4).Fix();
    fitter.Config().ParSettings(5).Fix();       
    fitter.Config().ParSettings(6).Fix();
    fitter.Config().ParSettings(7).Fix();

    fitter.Config().ParSettings(7).SetLimits(0,1);    
    
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit","Migrad");
    
// fiting FCN function directly (specify optionally data size and flag to indicate that is a chi2 fit)    
    fitter.FitFCN(8, globalChi2, 0, dataFIT_CHAN.Size()+dataFIT_REST.Size()+dataFIT_INT.Size()+dataFIT_NO_REACT.Size(), true);
    
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);


// plotting
    TCanvas *PROBABILITIES = new TCanvas("PROBABILITIES","PROBABILITIES", 0,0,900,900);
        
    PROBABILITIES->Divide(1,4);
    PROBABILITIES->cd(1);
        
        DATA_PROB_CHAN->SetTitle("state A");
        DATA_PROB_CHAN->GetXaxis()->SetTitle(" time (fs)");
        DATA_PROB_CHAN->GetXaxis()->CenterTitle();
        DATA_PROB_CHAN->GetYaxis()->SetTitle("Prob. (a.u.)");
        DATA_PROB_CHAN->GetYaxis()->CenterTitle();
        DATA_PROB_CHAN->Draw("apz");    
        
        FIT_CHAN->SetFitResult( result, iparCHAN);
        FIT_CHAN->SetRange(rangeFIT_CHAN().first, rangeFIT_CHAN().second);

        const int o_CHAN=2000;
        double Y_CHAN[o_CHAN], X_CHAN[o_CHAN];
        for(int h_CHAN=0; h_CHAN < o_CHAN; h_CHAN++)
        {
            X_CHAN[h_CHAN]=TIME.Min()+h_CHAN*(TIME.Max()-TIME.Min())/double(o_CHAN-1);
            
            Y_CHAN[h_CHAN]=FIT_CHAN->Eval(X_CHAN[h_CHAN]);
        
        }    
    
        TGraph *gr_CHAN=new TGraph(o_CHAN, X_CHAN, Y_CHAN);
        
        gr_CHAN->SetLineColor(4);
        gr_CHAN->SetLineWidth(3);
        gr_CHAN->Draw("L");           
//______________________________________________________________________________________________ 
    PROBABILITIES->cd(2);
        
        DATA_PROB_REST->SetTitle("state B");
        DATA_PROB_REST->GetXaxis()->SetTitle(" time (fs)");
        DATA_PROB_REST->GetXaxis()->CenterTitle();
        DATA_PROB_REST->GetYaxis()->SetTitle("Prob. (a.u.)");
        DATA_PROB_REST->GetYaxis()->CenterTitle();
        DATA_PROB_REST->Draw("apz");    
        
        FIT_REST->SetFitResult( result, iparREST);
        FIT_REST->SetRange(rangeFIT_REST().first, rangeFIT_REST().second);

        const int o_REST=2000;
        double Y_REST[o_REST], X_REST[o_REST];
        for(int h_REST=0; h_REST < o_REST; h_REST++)
        {
            X_REST[h_REST]=TIME.Min()+h_REST*(TIME.Max()-TIME.Min())/double(o_REST-1);
            
            Y_REST[h_REST]=FIT_REST->Eval(X_REST[h_REST]);
        
        }    
    
        TGraph *gr_REST=new TGraph(o_REST, X_REST, Y_REST);
        
        gr_REST->SetLineColor(4);
        gr_REST->SetLineWidth(3);
        gr_REST->Draw("L");      
//______________________________________________________________________________________________
        PROBABILITIES->cd(3);
        
        DATA_PROB_INT->SetTitle("state INT");
        DATA_PROB_INT->GetXaxis()->SetTitle(" time (fs)");
        DATA_PROB_INT->GetXaxis()->CenterTitle();
        DATA_PROB_INT->GetYaxis()->SetTitle("Prob. (a.u.)");
        DATA_PROB_INT->GetYaxis()->CenterTitle();
        DATA_PROB_INT->Draw("apz");    

        FIT_INT->SetFitResult( result, iparINT);
        FIT_INT->SetRange(rangeFIT_INT().first, rangeFIT_INT().second);

        const int o_INT=2000;
        double Y_INT[o_INT], X_INT[o_INT];
        for(int h_INT=0; h_INT<o_INT; h_INT++)
        {
            X_INT[h_INT]=TIME.Min()+h_INT*(TIME.Max()-TIME.Min())/double(o_INT-1);
            
            Y_INT[h_INT]=FIT_INT->Eval(X_INT[h_INT]);
        
        }    
    
        TGraph *gr_INT=new TGraph(o_INT, X_INT, Y_INT);
        
        gr_INT->SetLineColor(4);
        gr_INT->SetLineWidth(3);
        gr_INT->Draw("L");
//______________________________________________________________________________________________   
        PROBABILITIES->cd(4);
        
        DATA_PROB_NO_REACT->SetTitle("state PF");
        DATA_PROB_NO_REACT->GetXaxis()->SetTitle(" time (fs)");
        DATA_PROB_NO_REACT->GetXaxis()->CenterTitle();
        DATA_PROB_NO_REACT->GetYaxis()->SetTitle("Prob. (a.u.)");
        DATA_PROB_NO_REACT->GetYaxis()->CenterTitle();
        DATA_PROB_NO_REACT->Draw("apz");    

        FIT_NO_REACT->SetFitResult( result, iparNO_REACT);
        FIT_NO_REACT->SetRange(rangeFIT_NO_REACT().first, rangeFIT_NO_REACT().second);

        const int o_NO_REACT=2000;
        double Y_NO_REACT[o_NO_REACT], X_NO_REACT[o_NO_REACT];
        for(int h_NO_REACT=0; h_NO_REACT<o_NO_REACT; h_NO_REACT++)
        {
            X_NO_REACT[h_NO_REACT]=TIME.Min()+h_NO_REACT*(TIME.Max()-TIME.Min())/double(o_NO_REACT-1);
            
            Y_NO_REACT[h_NO_REACT]=FIT_NO_REACT->Eval(X_NO_REACT[h_NO_REACT]);
        
        }    
    
        TGraph *gr_NO_REACT=new TGraph(o_NO_REACT, X_NO_REACT, Y_NO_REACT);
        
        gr_NO_REACT->SetLineColor(4);
        gr_NO_REACT->SetLineWidth(3);
        gr_NO_REACT->Draw("L");
// ______________________________________________________________________________________________     
        const int o=2000;
        ofstream Results_lines("MULTIFIT-RESULTS_LINES");
        
        Results_lines << " TIME (fs) " << "\t"  << " state-SP " <<  "\t"  << " state-INT "  <<  "\t"  << " state-A " <<  "\t"  << " state-B " <<  endl;
        
        for(int h=0; h < o; h++)
        {
            
            Results_lines << X_CHAN[h] << "\t"  << Y_NO_REACT[h] <<  "\t"  << Y_INT[h]  <<  "\t"  << Y_CHAN[h] <<  "\t"  << Y_REST[h] <<  endl;
    
        }
        
        PROBABILITIES -> Print("4th-MultipleFit.pdf");

//  results and printing         
        ofstream Results("MULTIFIT-RESULTS_FIT");
        
        double a01=FIT_CHAN->GetParameter(0);
        double ea01=FIT_CHAN->GetParError(0);
        
        double b01=FIT_CHAN->GetParameter(1);
        double eb01=FIT_CHAN->GetParError(1);
        
        double c01=FIT_CHAN->GetParameter(2);
        double ec01=FIT_CHAN->GetParError(2);
        
        double d01=FIT_CHAN->GetParameter(3);
        double ed01=FIT_CHAN->GetParError(3);
        
        double e01=FIT_CHAN->GetParameter(4);
        double ee01=FIT_CHAN->GetParError(4);
        
        double f01=FIT_CHAN->GetParameter(5);
        double ef01=FIT_CHAN->GetParError(5);
        
        double tao=FIT_CHAN->GetParameter(6);
        double etao=FIT_CHAN->GetParError(6);
        
        double g01=FIT_CHAN->GetParameter(7);
        double eg01=FIT_CHAN->GetParError(7);
        
        double chi_red=FIT_CHAN->GetChisquare()/FIT_CHAN->GetNDF();
        
        double eChi_a01=ea01;
        double eChi_b01=eb01;
        double eChi_c01=ec01;
        double eChi_d01=ed01;
        double eChi_e01=ee01;
        double eChi_f01=ef01;
        double eChi_tao=etao;
        double eChi_g01=eg01;
        
        if(chi_red > 1)
            {
                eChi_a01=eChi_a01*sqrt(chi_red);
                eChi_b01=eChi_b01*sqrt(chi_red);
                eChi_c01=eChi_c01*sqrt(chi_red);
                eChi_d01=eChi_d01*sqrt(chi_red);
                eChi_e01=eChi_e01*sqrt(chi_red);
                eChi_f01=eChi_f01*sqrt(chi_red);
                eChi_tao=eChi_tao*sqrt(chi_red);
                eChi_g01=eChi_g01*sqrt(chi_red);
                        
            }
        
        
        Results.precision(6);
        
        Results << "\n chi_red  :\t " << chi_red  << "\t NDF  : \t"  << FIT_CHAN->GetNDF() <<  endl;
             
        Results << "\n kis         :\t " << a01 << "\t +/- \t"  << ea01 <<  "\t eChi \t"  << eChi_a01 <<  endl;
        Results << "\n ksi         :\t " << b01 << "\t +/- \t"  << eb01 << "\t eChi \t"  << eChi_b01 <<   endl;
        Results << "\n kas         :\t " << c01 << "\t +/- \t"  << ec01 << "\t eChi \t"  << eChi_c01 <<   endl;
        Results << "\n kai         :\t " << d01 << "\t +/- \t"  << ed01 << "\t eChi \t"  << eChi_d01 <<   endl;
        Results << "\n kas + kbs   :\t " << e01 << "\t +/- \t"  << ee01 << "\t eChi \t"  << eChi_e01 <<   endl;
        Results << "\n kai + kbi   :\t " << f01 << "\t +/- \t"  << ef01 << "\t eChi \t"  << eChi_f01 <<   endl;
        Results << "\n t0          :\t " << tao << "\t +/- \t"  << etao <<  "\t eChi \t"  << eChi_tao <<  endl;
        Results << "\n kba         :\t " << g01 << "\t +/- \t"  << eg01 <<  "\t eChi \t"  << eChi_g01 <<  endl;

        Results << "\n \t\t\t\t -------------------- \t \t " << endl;    


    
        timer.Stop();
    
       double cputime = timer.CpuTime();
       int hours=cputime/3600;
       int minutes=(cputime-3600*hours)/60;
       double seconds=cputime-3600*hours-60*minutes;
    cout << "\n\nCPU time (h:m:s)\t\t" << hours << ":" << minutes << ":" << seconds << "\n\n" << endl;
    
}



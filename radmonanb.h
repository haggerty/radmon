//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb  1 21:45:41 2013 by ROOT version 5.30/03
// from TTree r/radmon
// found on file: radmon.root
//////////////////////////////////////////////////////////

#ifndef radmonanb_h
#define radmonanb_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <limits>

#include <iostream>

class radmonanb {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           id;
   UInt_t          t;
   Int_t           channel;
   Float_t         i_n_set;
   Float_t         i_n;
   Float_t         v_n;
   Float_t         i_k_set;
   Float_t         i_k;
   Float_t         v_k;
   Float_t         i_s_set;
   Float_t         i_s;
   Float_t         v_s;
   Float_t         i_r_set;
   Float_t         i_r;
   Float_t         v_r;

   // List of branches
   TBranch        *b_id;   //!
   TBranch        *b_t;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_i_n_set;   //!
   TBranch        *b_i_n;   //!
   TBranch        *b_v_n;   //!
   TBranch        *b_i_k_set;   //!
   TBranch        *b_i_k;   //!
   TBranch        *b_v_k;   //!
   TBranch        *b_i_s_set;   //!
   TBranch        *b_i_s;   //!
   TBranch        *b_v_s;   //!
   TBranch        *b_i_r_set;   //!
   TBranch        *b_i_r;   //!
   TBranch        *b_v_r;   //!

   radmonanb(TTree *tree=0);
   virtual ~radmonanb();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop( Int_t idfirst = 0, Int_t idlast = std::numeric_limits<int>::max(), TString plotfile = "radmon.pdf" );
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef radmonanb_cxx
radmonanb::radmonanb(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("radmon.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("radmon.root");
      }
      f->GetObject("r",tree);

   }
   Init(tree);
}

radmonanb::~radmonanb()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t radmonanb::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t radmonanb::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void radmonanb::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("i_n_set", &i_n_set, &b_i_n_set);
   fChain->SetBranchAddress("i_n", &i_n, &b_i_n);
   fChain->SetBranchAddress("v_n", &v_n, &b_v_n);
   fChain->SetBranchAddress("i_k_set", &i_k_set, &b_i_k_set);
   fChain->SetBranchAddress("i_k", &i_k, &b_i_k);
   fChain->SetBranchAddress("v_k", &v_k, &b_v_k);
   fChain->SetBranchAddress("i_s_set", &i_s_set, &b_i_s_set);
   fChain->SetBranchAddress("i_s", &i_s, &b_i_s);
   fChain->SetBranchAddress("v_s", &v_s, &b_v_s);
   fChain->SetBranchAddress("i_r_set", &i_r_set, &b_i_r_set);
   fChain->SetBranchAddress("i_r", &i_r, &b_i_r);
   fChain->SetBranchAddress("v_r", &v_r, &b_v_r);
   Notify();
}

Bool_t radmonanb::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void radmonanb::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t radmonanb::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   std::cout << "entry: " << entry << std::endl;
   return 1;
}
#endif // #ifdef radmonanb_cxx

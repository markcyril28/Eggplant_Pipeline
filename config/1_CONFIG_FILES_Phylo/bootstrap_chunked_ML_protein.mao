; MEGA Analysis Options File
; Minimal ML Bootstrap Config (Chunk-compatible)
; Compatible with MEGA-CC v12

[ MEGAinfo ]
ver                                  = 12250409-x86_64 MS Windows          

[ DataSettings ]
datatype                             = snProtein
MissingBaseSymbol                    = ?
IdenticalBaseSymbol                  = .
GapSymbol                            = -

[ ProcessTypes ]
ppInfer                              = true
ppML                                 = true

[ AnalysisSettings ]
Analysis                             = Phylogeny Reconstruction
Statistical Method                   = Maximum Likelihood
Test of Phylogeny                    = Bootstrap
; Bootstrap Replicates omitted here, handled via CLI --bs_start / --bs_end
Substitutions Type                   = Amino acid
Model/Method                         = JTT with Freqs. (+F) model
Rates among Sites                    = Gamma Distributed With Invariant Sites (G+I)
No of Discrete Gamma Categories      = 5
Gaps/Missing Data                    = Use all sites
Tree Inference Options               = ====================
ML Heuristic Method                  = Subtree-Pruning-Regrafting - Extensive (SPR level 5)
Initial Tree for ML                  = NJ/MP (Default)
Branch Swap Filter                   = None
System Resource Usage                = ====================
Number of Threads                    = 1
Has Time Limit                       = False
Maximum Execution Time               = -1

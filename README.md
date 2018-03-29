# Subdistribution-beta-Stacy
The following code can be used to reproduced the analyses of the paper "Reinforced urns and the subdistribution beta-Stacy process prior for competing risks analysis" by Andrea Arfè, Stefano Peluso, and Pietro Muliere.

Code description:
Functions.R                        - contains the definitions of all functions required for the analyses
melanoma_gender_m1.R               - program for the analysis of Section 7 with m=1
melanoma_gender_m1000.R            - program for the analysis of Section 7 with m=1000
melanoma_gender_m1e6.R             - program for the analysis of Section 7 with m=1000000
melanoma_gender_m1_lognormal.R     - program for the analyses in the Supplementary Material, Appendix C, with m=1
melanoma_gender_m1000_lognormal.R  - program for the analyses in the Supplementary Material, Appendix C, with m=1
melanoma_gender_m1e6_lognormal.R   - program for the analyses in the Supplementary Material, Appendix C, with m=1000000
plot_pred_dist.R                   - plot of the predictive distributions for the analysis of Section 7
plot_pred_dist_lognormal.R         - plot of the predictive distributions for the analysis in the Supplementary Material, Appendix C
predictive_checks.R                - program to perform the predictive checks in the Supplementary Material, Appendix B
prior_simulations.R                - program to perform the analyses in the Supplementary Material, Appendix A
prior_standard_dev.R               - program to generate Figure 2 of Section 7

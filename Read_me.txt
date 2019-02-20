1) MAIN_gof_and_sparsity: 
"Starting from two MV-patterns (either simulated or real), this script can 
be used for estimating 1) the linear transformation by using the Tikhonov
regularization method, 2) the goodness-of-fit (GOF) and 3) the percentage of
sparsity (via Monte Carlo procedure that takes into account both the GOF
value and the rate of decay of the density curve (RDD)) as in Basti et
al. 2019. The results are then plotted, like in the PNG figure
called estimated_percentage_sparsity_original_value_75%."

2) MAIN_gof_and_deformation: 
"Starting from two MV-patterns (either simulated or real), this script can 
be used for estimating 1) the linear transformation by using the Tikhonov
regularization method, 2) the goodness-of-fit (GOF) and 3) the induced
pattern deformation (via Monte Carlo procedure which takes into account
both the GOF value and the rate of decay of the singular value curve(RDSV))
as in Basti et al. 2019. The results are then plotted, like in the PNG figure
called estimated_pattern_deformation_original_value_0.5."

3) The other functions are called in the those two scripts.

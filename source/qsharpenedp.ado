cap program drop qsharpenedp
program qsharpenedp
syntax varname
    * Sort the p-values in ascending order and generate a variable that codes
    * each p-value's rank
    quietly gen int original_sorting_order = _n
    quietly sort `1'
    quietly gen int rank = _n if `1'!=.

    quietly sum `1'
    local totalpvals = r(N)    

    * Set the initial counter to 1
    local qval = 1
    
    * Generate the variable that will contain the BKY (2006) sharpened q-values    
    qui: gen bky06_qval = 1 if `1'!=.
    
    * Set up a loop that begins by checking which hypotheses are rejected at
    * q = 1.000, then checks which hypotheses are rejected at q = 0.999, then
    * checks which hypotheses are rejected at q = 0.998, etc.  The loop ends by
    * checking which hypotheses are rejected at q = 0.001.
    
    while `qval' > 0 {
        * First Stage
        * Generate the adjusted first stage q level we are testing: q' = q/1+q
        local qval_adj = `qval'/(1+`qval')
        * Generate value q'*r/M
        qui: gen fdr_temp1 = `qval_adj'*rank/`totalpvals'
        * Generate binary variable checking condition p(r) <= q'*r/M
        qui: gen reject_temp1 = (fdr_temp1>=`1') if `1'!=.
        * Generate variable containing p-value ranks for all p-values that meet above condition
        qui: gen reject_rank1 = reject_temp1*rank
        * Record the rank of the largest p-value that meets above condition
        qui: egen total_rejected1 = max(reject_rank1)
    
        * Second Stage
        * Generate the second stage q level that accounts for hypotheses rejected 
        * in first stage: q_2st = q'*(M/m0)
        local qval_2st = `qval_adj'*(`totalpvals'/(`totalpvals'-total_rejected1[1]))
        * Generate value q_2st*r/M
        qui: gen fdr_temp2 = `qval_2st'*rank/`totalpvals'
        * Generate binary variable checking condition p(r) <= q_2st*r/M
        qui: gen reject_temp2 = (fdr_temp2>=`1') if `1'!=.
        * Generate variable containing p-value ranks for all p-values that meet above condition
        qui: gen reject_rank2 = reject_temp2*rank
        * Record the rank of the largest p-value that meets above condition
        qui: egen total_rejected2 = max(reject_rank2)
        
        * A p-value has been rejected at level q if its rank is less than or equal
        *to the rank of the max p-value that meets the above condition
        qui: replace bky06_qval = `qval' if rank <= total_rejected2 & rank~=.
        * Reduce q by 0.001 and repeat loop
        drop fdr_temp* reject_temp* reject_rank* total_rejected*
        local qval = `qval' - .001
    }
    quietly sort original_sorting_order
    display "Code has completed."
    display "Benjamini Krieger Yekutieli (2006) sharpened q-vals are in variable bky06_qval"
    display"Sorting order is the same as the original vector of p-values"
    drop original_sorting_order rank
end

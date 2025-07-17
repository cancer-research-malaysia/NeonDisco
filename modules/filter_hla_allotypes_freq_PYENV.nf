// 
process FILTER_HLA_ALLOTYPES_FREQ_PYENV {
    
    label 'filterHLAByFreq'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/5-PERCENT-TOP-COHORTWIDE-HLAS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        path(cohortWideHLAList)

    output:
        path("Frequency_filtered_HLA_allotypes.txt"), emit: freqFilteredHLAList
        path('HLA_allotypes_thresholded_cohortwide_string.txt'), emit: hlaAllotypesFreqFilteredString
        path("Frequency_filtering_report.txt")


    script:
    """
    echo "Running Python script to filter HLA types based on cohort-wide frequency..."
    if frequency-filter-HLAs--nf.py "${cohortWideHLAList}"; then
        echo "HLA allotype cohort-based frequency filtering completed successfully."
    else
        echo "Frequency filtering failed. Please check the logs for errors."
        exit 1
    fi
    """
    stub:
    """
    touch "Frequency_filtered_HLA_allotypes.txt"
    touch "HLA_allotypes_thresholded_cohortwide_string.txt"
    touch "Frequency_filtering_report.txt"
    echo "stub run finished!\thello my world!" > Frequency_filtered_HLA_allotypes.txt
    """
}

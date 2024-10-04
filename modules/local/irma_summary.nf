// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IRMA_SUMMARY {
  publishDir "${params.outdir}/", mode: params.publish_dir_mode
  label 'process_low'
  conda 'conda-forge::python=3.10 conda-forge::biopython=1.80 conda-forge::openpyxl=3.1.0 conda-forge::pandas=1.5.3 conda-forge::rich=12.6.0 conda-forge::typer=0.7.0 conda-forge::xlsxwriter=3.0.8 conda-forge::polars=0.17.9 conda-forge::pyarrow=11.0.0'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
    }

  input:
  path(irma_results)

  output:
  path('Ldt_Summary.txt'), emit: ldt_report
  path('Specimen_Stats.txt'), emit: specimen_report
  path('summary_fasta/**'), emit: fastas
  // path('summarize_irma_output.log'), emit: log

  script:
  """
  summarize_irma_output.py \\
    --coverage-depth-threshold $params.coverage_depth_threshold \\
    $irma_results
  # ln -s .command.log summarize_irma_output.log
  """
}

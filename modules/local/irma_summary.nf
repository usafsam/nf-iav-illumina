// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IRMA_SUMMARY {
  publishDir "${params.outdir}/", mode: params.publish_dir_mode
  label 'process_low'
  conda (params.enable_conda ? 'conda-forge::python=3.9 bioconda::biopython=1.78 conda-forge::openpyxl=3.0.7 conda-forge::pandas=1.2.4 conda-forge::rich=10.2.2 conda-forge::typer=0.3.2 conda-forge::xlsxwriter=1.4.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-693e24f156d01a5f55647120be99929b01b30949:609c862c3470382215fc1b2d9d8a4e9637b2e25f-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-693e24f156d01a5f55647120be99929b01b30949:609c862c3470382215fc1b2d9d8a4e9637b2e25f-0'
    }

  input:
  path(irma_results)

  output:
  path('Ldt_Summary.txt'), emit: ldt_report
  path('Specimen_Stats.txt'), emit: specimen_report
  path('summary_fasta/**'), emit: fastas
  path('summarize_irma_output.log'), emit: log

  script:
  """
  summarize_irma_output.py \\
    --coverage-depth-threshold $params.coverage_depth_threshold \\
    $irma_results
  ln -s .command.log summarize_irma_output.log
  """
}

$proteinData = Import-Csv -Path 'C:\Users\jtzve\Desktop\Sufo_Tyrosine\04_MotifEnrichment\out\protein_ids_for_FASTA_download.tsv' -Delimiter "`t"
foreach ($row in $proteinData) {
    $id = $row.ID # The column name is 'ID' as per your description
    $url = "https://www.uniprot.org/uniprot/$id.fasta"
    $outputFile = "C:\Users\jtzve\Desktop\Sufo_Tyrosine\04_MotifEnrichment\data\SwissProt_FASA_files\$id.fasta"
    Invoke-WebRequest -Uri $url -OutFile $outputFile
}

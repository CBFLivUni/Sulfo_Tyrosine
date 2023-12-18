$source = 'C:\Users\jtzve\Desktop\Sufo_Tyrosine\02_Sulfotyrosine_peptidoform_aggreagation'
$destination = 'C:\Users\jtzve\Desktop\Sufo_Tyrosine\01_Sulfotyrosine_Data_Processing\data_to_process'

# Create the destination directory if it doesn't exist
if (-not (Test-Path -Path $destination)) {
    New-Item -ItemType Directory -Path $destination
}

# move the files and folders, only keeping .xml
Get-ChildItem -Path $source -Recurse | Where-Object { -not $_.PSIsContainer -and $_.Extension -eq '.xml' } | ForEach-Object {
    $dest = $_.FullName.Replace($source, $destination)
    if (-not (Test-Path -Path (Split-Path -Path $dest -Parent))) {
        New-Item -ItemType Directory -Path (Split-Path -Path $dest -Parent)
    }
    Move-Item -Path $_.FullName -Destination $dest
}

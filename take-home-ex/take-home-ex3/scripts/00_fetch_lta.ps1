# =============================================================================
# 00_fetch_lta.ps1 -- monthly passenger volumes from LTA DataMall (Windows)
#
# The PowerShell twin of 00_fetch_lta.sh. PV/ODTrain gives trips between pairs of
# stations by day type and hour; PV/Train gives tap-ins and tap-outs by station,
# and is used only to check the totals built from the OD file. DataMall keeps the
# three latest months, and each request returns a link that expires in 5 minutes.
#
# The key is asked for at a hidden prompt, so it never lands in the repository or
# in the shell history. Files are saved to data/raw/lta/.
# =============================================================================
param([string[]]$Months = @("202607", "202608", "202609"))

$root = Split-Path $PSScriptRoot -Parent
$out  = Join-Path $root "data\raw\lta"
New-Item -ItemType Directory -Force -Path $out | Out-Null

$key = $env:LTA_KEY
if (-not $key) {
    $secure = Read-Host "DataMall AccountKey" -AsSecureString
    $bstr   = [Runtime.InteropServices.Marshal]::SecureStringToBSTR($secure)
    $key    = [Runtime.InteropServices.Marshal]::PtrToStringAuto($bstr)
    [Runtime.InteropServices.Marshal]::ZeroFreeBSTR($bstr)
}
if (-not $key) { Write-Host "No key given."; exit 1 }

[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$ProgressPreference = "SilentlyContinue"
$base = "https://datamall2.mytransport.sg/ltaodataservice"
$headers = @{ AccountKey = $key; accept = "application/json" }

foreach ($month in $Months) {
    foreach ($endpoint in @("PV/ODTrain", "PV/Train")) {
        $name = ($endpoint -replace "/", "") + "_" + $month + ".zip"
        $url  = $base + "/" + $endpoint + "?Date=" + $month
        try {
            $response = Invoke-RestMethod -Uri $url -Headers $headers -Method Get
            $link = $response.value[0].Link
            if (-not $link) {
                Write-Host ("FAILED " + $endpoint + " " + $month + ": no link in response: " + ($response | ConvertTo-Json -Compress))
                continue
            }
            $file = Join-Path $out $name
            Invoke-WebRequest -Uri $link -OutFile $file -UseBasicParsing
            Write-Host ("saved " + $name + " (" + (Get-Item $file).Length + " bytes)")
        } catch {
            Write-Host ("FAILED " + $endpoint + " " + $month + ": " + $_.Exception.Message)
        }
    }
}

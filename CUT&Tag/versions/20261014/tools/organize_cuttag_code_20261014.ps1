param(
    [string]$Workspace = 'C:\Users\Luca\Documents\Codex\2026-09-17\cut-tag',
    [string]$PackageName = 'CUTTag_code_versions_20261014'
)

$ErrorActionPreference = 'Stop'
$workspaceFull = [IO.Path]::GetFullPath($Workspace).TrimEnd([IO.Path]::DirectorySeparatorChar)
if (-not (Test-Path -LiteralPath $workspaceFull -PathType Container)) { throw 'Workspace not found' }
if ($PackageName -notmatch '^CUTTag_code_versions_[0-9]{8}$') { throw 'Unexpected package name' }
$outputsFull = Join-Path $workspaceFull 'outputs'
$packageFull = [IO.Path]::GetFullPath((Join-Path $outputsFull $PackageName))
if (-not $packageFull.StartsWith($outputsFull + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) {
    throw 'Package is outside the intended outputs directory'
}
if (-not (Test-Path -LiteralPath (Join-Path $packageFull 'README.md'))) {
    throw 'Create the authored README and docs with apply_patch before building'
}
$utf8 = [Text.UTF8Encoding]::new($false)

function GeneratedText([string]$Path, [string]$Text) {
    # Generated indexes/metadata only; never change analysis source files.
    [IO.Directory]::CreateDirectory([IO.Path]::GetDirectoryName($Path)) | Out-Null
    [IO.File]::WriteAllText($Path, $Text.Replace("`r`n", "`n").Replace("`r", "`n"), $utf8)
}
function ArchiveCopy([string]$Source, [string]$Destination) {
    if (-not (Test-Path -LiteralPath $Source -PathType Leaf)) { throw "Missing source: $Source" }
    [IO.Directory]::CreateDirectory([IO.Path]::GetDirectoryName($Destination)) | Out-Null
    $sourceHash = (Get-FileHash -LiteralPath $Source -Algorithm SHA256).Hash.ToLowerInvariant()
    if (Test-Path -LiteralPath $Destination) {
        if ((Get-FileHash -LiteralPath $Destination -Algorithm SHA256).Hash.ToLowerInvariant() -ne $sourceHash) {
            throw "Existing archived file differs; refusing to overwrite: $Destination"
        }
    } else {
        Copy-Item -LiteralPath $Source -Destination $Destination
    }
    $destinationHash = (Get-FileHash -LiteralPath $Destination -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($sourceHash -ne $destinationHash) { throw "Hash mismatch: $Source" }
    return $sourceHash
}

$groups = [ordered]@{
    '00_environment' = @{ Title='环境准备与离线预热'; Note='prewarm 默认检查字面量 Conda 规格、conda-meta 和 SEACR 执行文件。安装模式应在有网管理节点运行。Nextflow 内部 hash 依赖所用版本。缺失的 cutandtag_env.sh/local_conda.config 见总说明。' }
    '01_v1_legacy' = @{ Title='旧结果 v1 的轨道补充'; Note='历史脚本读取 results/，gsize=2.7e9，目标为单重复和 pooled 共 12 条 MACS2 ppois 轨道。旧完整主流程提交文件未找到。以 v2 为主时不要把这个脚本当作 v2 后处理。' }
    '02_v2_main' = @{ Title='猪参考主流程 v2、轨道与统计'; Note='run→postprocess→TSS 是输入依赖顺序；collect 在需要统计的结果生成后运行。只有 collect 没有 --check 模式。test_collect 与 collect 保持同目录。参数、缓存和重复解释见 docs。' }
    '03_v2_hypothesis_screen' = @{ Title='严格现有峰的候选筛查'; Note='不放宽阈值、不重做 BAM。生成 peak_sets、candidate/rRNA/5S/tRNA 注释和共占位，供后续定量与敏感性脚本使用。候选基因邻近不等于功能调控或蛋白互作。' }
    '04_signal_v1_legacy' = @{ Title='初版区域信号定量与审计'; Note='multiBigwigSummary 版本保留供历史溯源。区域均值的缺失碱基/分母解释在精确定量 v2 修正。本地审计代码依赖原工作区相对结构，移动后需调整 ROOT/DATA/SCREEN/OUT。' }
    '05_signal_v2_exact' = @{ Title='精确信号定量 v2'; Note='Python 和 Slurm 必须同时部署。wrapper 调用 $PROJECT/config/quantify_cuttag_signal_qc_v2.py；测试与 Python 放同目录。该目录与旧交付包中的同名 Python/Slurm 是相同 SHA256，不额外制造新版本。历史方法说明中的尚未执行等状态，以本归档版本表更新。' }
    '06_sensitivity_q05' = @{ Title='narrow q0.05 参数敏感性'; Note='job 51897 对应脚本；Tag 两重复分别调用，保留同状态 IgG 和局部背景。共同区间原则不变，不追加 q0.01 二次过滤，不生成新 bigWig。' }
    '07_sensitivity_broad_nolambda' = @{ Title='broad 与 nolambda 参数敏感性'; Note='同一脚本用 broad/nolambda 参数分别提交。broad 输出 broadPeak/gappedPeak，nolambda 输出 narrowPeak。比较集固定为 q05 job-51897；job 52126 是已提供的 nolambda 结果，broad 完成结果尚未提供。' }
    '08_local_tables' = @{ Title='本地统计表工具的版本'; Note='v1_template 是模板生成；v1_filled 是填表与解析/检查；fragment_draft 使用 tmp/预置资料，为历史草稿；fragment_preserve 是读取完成 TSV 并保留公式缓存的后续版；v2 为新版参数说明与统计填数。inspection 是早期参考表查看工具。多数脚本写死路径，依赖 @oai/artifact-tool 或 openpyxl，运行前读代码中的输入输出。' }
    '09_local_reports' = @{ Title='2026-10-13 本地报告生成代码'; Note='v1 生成简洁汇报；v2_focus 使用模板并突出区域/基因结果，执行顺序是 prepare_data→inspect_source→build_focus。本归档只保存代码，未携带模板、输入数据、图片和 artifact runtime。' }
}

$entries = [Collections.Generic.List[object]]::new()
function AddEntry([string]$Source, [string]$Target, [string]$Purpose, [string]$Status) {
    $entries.Add([pscustomobject]@{ Source=$Source; Target=$Target; Purpose=$Purpose; Status=$Status })
}
AddEntry 'prewarm_cutandtag_envs.nf' '00_environment/prewarm_cutandtag_envs.nf' '离线环境预热/检查' '历史修正后版本'
AddEntry 'outputs/macs2_igg_pval_tracks.slurm' '01_v1_legacy/macs2_igg_pval_tracks.slurm' '旧结果单重复+pooled ppois 轨道' '历史；成功全程未确认'
AddEntry 'outputs/run_cutandrun_pig_v2.slurm' '02_v2_main/run_cutandrun_pig_v2.slurm' '猪参考 v2 nf-core driver' '44536 成功'
AddEntry 'outputs/postprocess_cutandtag_pig_v2.slurm' '02_v2_main/postprocess_cutandtag_pig_v2.slurm' '四条 pooled ppois 与两 caller 共同峰' '44538 对应后处理'
AddEntry 'outputs/run_cuttag_tss_v2.slurm' '02_v2_main/run_cuttag_tss_v2.slurm' '全基因 TSS coverage profiles' '脚本存在；完成状态待日志核对'
AddEntry 'outputs/collect_cuttag_v2_stats.slurm' '02_v2_main/collect_cuttag_v2_stats.slurm' '收集 v2 统计至 JSON.GZ' '44539 对应任务'
AddEntry 'outputs/test_collect_cuttag_v2_stats.py' '02_v2_main/test_collect_cuttag_v2_stats.py' '统计收集器既有测试' '与 Slurm 同目录'
AddEntry 'outputs/screen_cuttag_interaction_hypotheses.slurm' '03_v2_hypothesis_screen/screen_cuttag_interaction_hypotheses.slurm' '严格峰注释和共占位筛查' '已有结果包'
AddEntry 'outputs/quantify_cuttag_hypothesis_signals.slurm' '04_signal_v1_legacy/quantify_cuttag_hypothesis_signals.slurm' '初版 multiBigwigSummary 区域定量' '被精确定量 v2 替代'
AddEntry 'work/audit_hypothesis_signal_quantification.py' '04_signal_v1_legacy/local_audit/audit_hypothesis_signal_quantification.py' '本地区域信号结果审计' '相对路径依赖原工作区'
AddEntry 'outputs/quantify_cuttag_signal_qc_v2.py' '05_signal_v2_exact/quantify_cuttag_signal_qc_v2.py' '原始 interval 积分与匹配比较' '51892 对应代码'
AddEntry 'outputs/quantify_cuttag_signal_qc_v2.slurm' '05_signal_v2_exact/quantify_cuttag_signal_qc_v2.slurm' '精确定量 v2 Slurm wrapper' '51892 对应代码'
AddEntry 'outputs/test_quantify_cuttag_signal_qc_v2.py' '05_signal_v2_exact/test_quantify_cuttag_signal_qc_v2.py' '精确定量既有回归测试' '与 helper 同目录'
AddEntry 'outputs/cuttag_macs2_q05_sensitivity.slurm' '06_sensitivity_q05/cuttag_macs2_q05_sensitivity.slurm' 'Tag narrow q0.05 与严格集合对照' '51897 成功'
AddEntry 'outputs/cuttag_macs2_broad_nolambda_sensitivity.slurm' '07_sensitivity_broad_nolambda/cuttag_macs2_broad_nolambda_sensitivity.slurm' 'broad/nolambda 两模式独立分析' 'nolambda 52126 成功；broad 未核验'
AddEntry 'work/build_cuttag_macs2_stats.mjs' '08_local_tables/v1_template/build_cuttag_macs2_stats.mjs' 'MACS2 项目统计模板生成' '本地历史模板'
AddEntry 'work/fill_workbook1_cuttag.mjs' '08_local_tables/v1_filled/fill_workbook1_cuttag.mjs' '原 Sheet3 填统计' '本地已填表对应版本'
AddEntry 'work/parse_cuttag_inputs.py' '08_local_tables/v1_filled/parse_cuttag_inputs.py' '解析旧统计输入文本' '本地辅助'
AddEntry 'work/inspect_multimqc_blocks.py' '08_local_tables/v1_filled/inspect_multimqc_blocks.py' '检查 MultiQC 统计字段' '本地辅助'
AddEntry 'work/inspect_current_workbook.py' '08_local_tables/v1_filled/inspect_current_workbook.py' '检查工作簿及数值' '本地辅助'
AddEntry 'work/inspect_workbook1_sheet3.mjs' '08_local_tables/v1_filled/inspect_workbook1_sheet3.mjs' 'Sheet3 片段/IgG 列预览' '本地辅助'
AddEntry 'work/update_fragment_igg.mjs' '08_local_tables/v1_fragment_draft/update_fragment_igg.mjs' 'tmp 与预置值的片段/IgG 填表草稿' '历史草稿；用 preserve 版溯源正式填表'
AddEntry 'work/update_fragment_igg_preserve.py' '08_local_tables/v1_fragment_preserve/update_fragment_igg_preserve.py' '完成 TSV 的片段均值/IgG 命名，保留缓存' '后续修订'
AddEntry 'work/check_fragment_igg_workbook.py' '08_local_tables/v1_fragment_preserve/check_fragment_igg_workbook.py' '核对片段/IgG 填表前后' '本地辅助'
AddEntry 'outputs/update_cuttag_v2_workbook.mjs' '08_local_tables/v2/update_cuttag_v2_workbook.mjs' 'v2 参数说明、--fill-v2 新统计' '依赖历史 JSON.GZ 与工作簿'
AddEntry 'work/inspect_reference.mjs' '08_local_tables/inspection/inspect_reference.mjs' '早期参考工作簿查看' '历史辅助；只保存代码'
AddEntry 'work/artifact/inspect_reference.mjs' '08_local_tables/inspection/artifact/inspect_reference.mjs' '早期 artifact 参考工作簿查看' '另一个内容不同的原版本'
AddEntry 'work/ppt_report_20261013/build_report.mjs' '09_local_reports/20261013_v1/build_report.mjs' '简洁汇报 PPT 生成' '本地报告版本'
AddEntry 'work/ppt_focus_20261013/prepare_data.mjs' '09_local_reports/20261013_v2_focus/prepare_data.mjs' '重点报告数据准备' '本地报告版本'
AddEntry 'work/ppt_focus_20261013/inspect_source.mjs' '09_local_reports/20261013_v2_focus/inspect_source.mjs' '模板/旧汇报检查' '本地报告版本'
AddEntry 'work/ppt_focus_20261013/build_focus.mjs' '09_local_reports/20261013_v2_focus/build_focus.mjs' '模板版区域/基因报告生成' '本地报告版本'
AddEntry 'outputs/CUTTag背景核查与敏感性分析说明.md' '05_signal_v2_exact/历史方法说明.md' '保留当时信号 QC 与敏感性说明' '历史状态文字以总版本表更新'
AddEntry 'outputs/hypothesis_signal_audit_20261009/结果解读与定量修正说明.md' '04_signal_v1_legacy/local_audit/历史审计说明.md' '初版定量修正记录' '历史方法记录'

$manifestEntries = [Collections.Generic.List[object]]::new()
foreach ($entry in $entries) {
    $src = Join-Path $workspaceFull $entry.Source
    $dst = Join-Path $packageFull $entry.Target
    $sha = ArchiveCopy $src $dst
    $item = Get-Item -LiteralPath $src
    $newlines = $null
    if ($entry.Target -match '\.(slurm|sh|nf|py|mjs|ps1)$') {
        $content = [IO.File]::ReadAllText($src)
        $crlf = ([regex]::Matches($content, "`r`n")).Count
        $lf = ([regex]::Matches($content, "(?<!`r)`n")).Count
        $newlines = if ($crlf -and $lf) { 'mixed' } elseif ($crlf) { 'CRLF' } elseif ($lf) { 'LF' } else { 'no_newline' }
    }
    $manifestEntries.Add([pscustomobject]@{
        source=$entry.Source; archived=$entry.Target; purpose=$entry.Purpose; status=$entry.Status
        bytes=$item.Length; sha256=$sha; newline_style=$newlines; copy_matches_source=$true
    })
}

# Preserve the local third-party snapshot, including LICENSE and dotfiles.
$upstreamSource = Join-Path $workspaceFull 'work/cutandrun-source/cutandrun-3.2.2'
$upstreamFiles = @(Get-ChildItem -LiteralPath $upstreamSource -File -Force -Recurse)
foreach ($file in $upstreamFiles) {
    $rel = [IO.Path]::GetRelativePath($upstreamSource, $file.FullName).Replace('\','/')
    if ($rel.StartsWith('.git/')) { continue }
    $target = '90_nfcore_source/cutandrun-3.2.2/' + $rel
    $sha = ArchiveCopy $file.FullName (Join-Path $packageFull $target)
    $manifestEntries.Add([pscustomobject]@{
        source='work/cutandrun-source/cutandrun-3.2.2/' + $rel; archived=$target
        purpose='本地 nf-core/cutandrun 3.2.2 源码快照'; status='未与服务器实际源码逐文件比较'
        bytes=$file.Length; sha256=$sha; newline_style=$null; copy_matches_source=$true
    })
}
GeneratedText (Join-Path $packageFull '90_nfcore_source/README.md') "# 本地 nf-core/cutandrun 3.2.2 源码快照`n`n保留本地 work/cutandrun-source/cutandrun-3.2.2 的全部 $($upstreamFiles.Count) 个文件，包括 LICENSE、modules、subworkflows、schema 和 docs。快照与服务器实际 pipeline 是否完全一致尚未逐文件比较。不要通过覆盖服务器 pipeline 来完成归档。`n"

$builderTarget='tools/organize_cuttag_code_20261014.ps1'
$builderHash=ArchiveCopy $PSCommandPath (Join-Path $packageFull $builderTarget)
$manifestEntries.Add([pscustomobject]@{
    source=[IO.Path]::GetRelativePath($workspaceFull,$PSCommandPath).Replace('\','/'); archived=$builderTarget
    purpose='代码归档、来源清单和 SHA256 生成'; status='本次打包工具；需要原工作区'
    bytes=(Get-Item -LiteralPath $PSCommandPath).Length; sha256=$builderHash; newline_style='LF'; copy_matches_source=$true
})

$inventory = [Text.StringBuilder]::new()
[void]$inventory.AppendLine('# 自定义代码清单')
[void]$inventory.AppendLine('')
[void]$inventory.AppendLine('代码保持原文件名，按目录分类；所有来源和 SHA256 见 manifest.json。旧交付包中重复的精确定量 Python/Slurm 与本包相同文件去重，记录在 aliases。')
foreach ($key in $groups.Keys) {
    $group = $groups[$key]
    $matches = @($entries | Where-Object { $_.Target.StartsWith($key+'/') })
    $readme = [Text.StringBuilder]::new()
    [void]$readme.AppendLine('# ' + $group.Title)
    [void]$readme.AppendLine('')
    [void]$readme.AppendLine($group.Note)
    [void]$readme.AppendLine('')
    [void]$readme.AppendLine('总运行说明：[运行顺序与输入输出](../docs/运行顺序与输入输出.md)。')
    [void]$readme.AppendLine('')
    [void]$readme.AppendLine('|文件|用途|状态|')
    [void]$readme.AppendLine('|---|---|---|')
    [void]$inventory.AppendLine('')
    [void]$inventory.AppendLine('## ' + $key + '：' + $group.Title)
    [void]$inventory.AppendLine('')
    [void]$inventory.AppendLine('|归档文件|用途|原始位置|状态|')
    [void]$inventory.AppendLine('|---|---|---|---|')
    foreach ($entry in $matches) {
        $sub=$entry.Target.Substring($key.Length+1)
        [void]$readme.AppendLine('|['+$sub+']('+ $sub +')|'+$entry.Purpose+'|'+$entry.Status+'|')
        [void]$inventory.AppendLine('|['+$entry.Target+']('+ $entry.Target +')|'+$entry.Purpose+'|`'+$entry.Source+'`|'+$entry.Status+'|')
    }
    GeneratedText (Join-Path $packageFull "$key/README.md") $readme.ToString()
}
[void]$inventory.AppendLine('')
[void]$inventory.AppendLine('上游源码单独存放在 90_nfcore_source/cutandrun-3.2.2；自定义归档工具在 tools/。不含结果数据或软件运行环境。')
GeneratedText (Join-Path $packageFull '文件清单.md') $inventory.ToString()

$aliases=@()
foreach($pair in @(
    @('outputs/cuttag_signal_qc_v2_20261009/quantify_cuttag_signal_qc_v2.py','outputs/quantify_cuttag_signal_qc_v2.py'),
    @('outputs/cuttag_signal_qc_v2_20261009/quantify_cuttag_signal_qc_v2.slurm','outputs/quantify_cuttag_signal_qc_v2.slurm'),
    @('outputs/cuttag_signal_qc_v2_20261009/CUTTag背景核查与敏感性分析说明.md','outputs/CUTTag背景核查与敏感性分析说明.md')
)) {
    $a=Join-Path $workspaceFull $pair[0];$b=Join-Path $workspaceFull $pair[1]
    if (-not (Test-Path -LiteralPath $a)) { continue }
    $same=(Get-FileHash -LiteralPath $a).Hash -eq (Get-FileHash -LiteralPath $b).Hash
    if (-not $same) { throw "Alias content differs: $a" }
    $aliases += @{ source=$pair[0]; same_as=$pair[1]; sha256_matches=$same }
}
$customScripts=@($entries|Where-Object {$_.Target -match '\.(slurm|sh|nf|py|mjs)$'})
$manifest=[ordered]@{
    package=$PackageName; archive_date='2026-10-14'; timezone='Asia/Shanghai'; source_workspace=$workspaceFull
    policy='Copies preserve source bytes; existing analysis files are not moved, deleted, renamed, or edited.'
    custom_analysis_and_local_scripts=$customScripts.Count; packaging_tools=1
    upstream_snapshot_files=$upstreamFiles.Count; aliases=$aliases; files=$manifestEntries.ToArray()
    missing_server_files=@('v1 run_cutandrun.slurm','config/params.yaml','config/samplesheet.csv','cutandtag_env.sh','local_conda.config','fragment length Slurm source','pre-v2 merged bigWig script')
}
GeneratedText (Join-Path $packageFull 'manifest.json') ($manifest | ConvertTo-Json -Depth 8)

# An inventory check detects accidentally omitted locally available custom code.
$customPaths=@(Get-ChildItem -LiteralPath $workspaceFull -File | Where-Object {$_.Extension -in '.nf','.sh','.slurm','.py','.mjs','.ps1','.config','.yaml','.yml','.R'})
$customPaths+=@(Get-ChildItem -LiteralPath (Join-Path $workspaceFull 'outputs') -File | Where-Object {$_.Extension -in '.nf','.sh','.slurm','.py','.mjs','.ps1','.config','.yaml','.yml','.R'})
$customPaths+=@(Get-ChildItem -LiteralPath (Join-Path $workspaceFull 'work') -File | Where-Object {$_.Extension -in '.nf','.sh','.slurm','.py','.mjs','.ps1','.config','.yaml','.yml','.R'})
foreach($d in 'work/artifact','work/ppt_report_20261013','work/ppt_focus_20261013') {
    $customPaths+=@(Get-ChildItem -LiteralPath (Join-Path $workspaceFull $d) -File | Where-Object {$_.Extension -in '.py','.mjs','.ps1'})
}
$recorded=@($manifestEntries|ForEach-Object{$_.source})
$omitted=@($customPaths | ForEach-Object { [IO.Path]::GetRelativePath($workspaceFull,$_.FullName).Replace('\','/') } | Where-Object {$_ -notin $recorded})
if($omitted.Count) { throw ('Custom code omitted: '+($omitted -join ', ')) }
Write-Output "Packaged $($customScripts.Count) custom scripts + 1 packaging tool; $($upstreamFiles.Count) upstream files. Source hashes match."
Write-Output $packageFull

$validation=[ordered]@{
    package=$PackageName; validated_at=[DateTimeOffset]::Now.ToString('o')
    source_copy_sha256='passed'; custom_code_inventory='passed'; upstream_license_present=(Test-Path -LiteralPath (Join-Path $packageFull '90_nfcore_source/cutandrun-3.2.2/LICENSE'))
    syntax_checks=@(); python_regression_tests='not rerun during packaging; original test sources are included'
    nextflow_compile='not run; cluster Nextflow and Conda cache are external dependencies'
    server_execution='not performed'
    zip_integrity_report='../'+$PackageName+'.zip.validation.json'
}
$nodeCmd=Get-Command node -ErrorAction SilentlyContinue
$bashCandidates=@('D:\游览器下载\Git\bin\bash.exe','C:\Program Files\Git\bin\bash.exe')
$bashExe=$bashCandidates | Where-Object {Test-Path -LiteralPath $_ -PathType Leaf} | Select-Object -First 1
foreach($entry in $customScripts) {
    $file=Join-Path $packageFull $entry.Target
    if($entry.Target.EndsWith('.mjs')) {
        if($nodeCmd) {
            $messages=@(& $nodeCmd.Source --check $file 2>&1)
            if($LASTEXITCODE -ne 0){throw ('Node syntax failed: '+$entry.Target+' '+($messages -join ' '))}
            $validation.syntax_checks+=@{file=$entry.Target; tool='node --check'; status='passed'}
        } else {$validation.syntax_checks+=@{file=$entry.Target; tool='node --check'; status='skipped: node unavailable'}}
    } elseif($entry.Target.EndsWith('.slurm')) {
        if($bashExe) {
            $messages=@(& $bashExe -n $file.Replace('\','/') 2>&1)
            if($LASTEXITCODE -ne 0) {
                $failure=$messages -join ' '
                if($failure -match 'NtCreateDirectoryObject') {
                    # MSYS cannot initialize in this sandbox; this is not a script syntax result.
                    $bashExe=$null
                    $validation.syntax_checks+=@{file=$entry.Target; tool='bash -n'; status='skipped: Git Bash runtime blocked by local sandbox'}
                } else {throw ('Bash syntax failed: '+$entry.Target+' '+$failure)}
            } else {$validation.syntax_checks+=@{file=$entry.Target; tool='bash -n'; status='passed'}}
        } else {$validation.syntax_checks+=@{file=$entry.Target; tool='bash -n'; status='skipped: Git Bash unavailable or runtime blocked'}}
    }
}
GeneratedText (Join-Path $packageFull 'validation.json') ($validation | ConvertTo-Json -Depth 7)

$sumLines=[Collections.Generic.List[string]]::new()
$packageFiles=@(Get-ChildItem -LiteralPath $packageFull -Recurse -File -Force | Where-Object {$_.Name -ne 'SHA256SUMS.txt'} | Sort-Object FullName)
foreach($file in $packageFiles) {
    $relative=[IO.Path]::GetRelativePath($packageFull,$file.FullName).Replace('\','/')
    $sha=(Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    $sumLines.Add($sha+'  '+$relative)
}
GeneratedText (Join-Path $packageFull 'SHA256SUMS.txt') (($sumLines -join "`n")+"`n")

Add-Type -AssemblyName System.IO.Compression
Add-Type -AssemblyName System.IO.Compression.FileSystem
$archiveFull=Join-Path $outputsFull ($PackageName+'.zip')
$archivePartial=$archiveFull+'.partial'
if((Test-Path -LiteralPath $archiveFull) -or (Test-Path -LiteralPath $archivePartial)) {
    throw 'Archive or partial archive already exists; refusing overwrite'
}
[IO.Compression.ZipFile]::CreateFromDirectory($packageFull,$archivePartial,[IO.Compression.CompressionLevel]::Optimal,$true)
$expected=@{}
foreach($file in (Get-ChildItem -LiteralPath $packageFull -Recurse -File -Force)) {
    $relative=$PackageName+'/'+[IO.Path]::GetRelativePath($packageFull,$file.FullName).Replace('\','/')
    $expected[$relative]=(Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
}
$zip=[IO.Compression.ZipFile]::OpenRead($archivePartial)
$verified=0
try {
    $seen=@{}
    foreach($entry in $zip.Entries) {
        if($entry.FullName.EndsWith('/')){continue}
        if(-not $expected.ContainsKey($entry.FullName) -or $seen.ContainsKey($entry.FullName)){throw ('Unexpected/duplicate ZIP entry: '+$entry.FullName)}
        $stream=$entry.Open()
        $hasher=[Security.Cryptography.SHA256]::Create()
        try {$actual=[Convert]::ToHexString($hasher.ComputeHash($stream)).ToLowerInvariant()}
        finally {$hasher.Dispose();$stream.Dispose()}
        if($actual -ne $expected[$entry.FullName]){throw ('ZIP content SHA256 mismatch: '+$entry.FullName)}
        $seen[$entry.FullName]=$true
        $verified++
    }
    if($verified -ne $expected.Count){throw 'ZIP has missing files'}
} finally {$zip.Dispose()}
Move-Item -LiteralPath $archivePartial -Destination $archiveFull
$archiveHash=(Get-FileHash -LiteralPath $archiveFull -Algorithm SHA256).Hash.ToLowerInvariant()
GeneratedText ($archiveFull+'.sha256.txt') ($archiveHash+'  '+[IO.Path]::GetFileName($archiveFull)+"`n")
GeneratedText ($archiveFull+'.validation.json') (@{
    package=$PackageName; zip_file_count=$verified; all_zip_entry_sha256_match=$true
    archive_sha256=$archiveHash; bytes=(Get-Item -LiteralPath $archiveFull).Length
    validated_at=[DateTimeOffset]::Now.ToString('o'); original_analysis_files_modified=$false
} | ConvertTo-Json)
Write-Output "Syntax checks: $($validation.syntax_checks.Count); ZIP entries verified: $verified"
Write-Output $archiveFull

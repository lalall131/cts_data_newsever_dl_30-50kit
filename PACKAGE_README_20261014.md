# CTS-Atlas 8类Kit下载与输入准备包（2026-10-14）

本版在原来的5类人类批次上，新增3类小鼠批次。完整名单为 **278个Source data ID/GSM测序样本、339个Run**。原5类名单保留，Source data ID/GSM是测序记录计数单位，不代表每个记录都是独立生物学重复，也不是单细胞的细胞数。

## 本次新增下载

| Kit | 物种 | GEO | 入选GSM | Run | 原始归档FASTQ.gz文件 | ENA已知压缩体积 |
|---|---|---|---:|---:|---:|---:|
| SPI-PBLR | Mus musculus | GSE144205 | 14 | 14 | 14 | 0.831 GB |
| scFLUENT-seq SPI | Mus musculus | GSE278777 | 12 | 12 | 24 | 548.058 GB |
| scGRO-seq SR | Mus musculus | GSE242176 | 50 | 50 | 82 | 11.967 GB |
| 合计 | 小鼠 | 3项研究 | **76** | **76** | **120** | **560.856 GB** |

这里是ENA列出的原始归档gzip体积，下载使用这些原文件，因而可以比对官方大小和MD5。不是SRA缓存体积、解压后体积、未来BAM体积或最终项目磁盘占用。另有26个未入选的合适scGRO候选，名单也保留。

### SPI-PBLR

GSE144205的PacBio RS II/Sequel小鼠MEL染色质RNA长读数据，全部14个测序样本均保留。

- 8个全基因组LRS样本：未诱导/红系分化诱导、各两个重复、RSII和Sequel记录。4对相同条件/重复编号的不同平台记录单独保留；不能当作8个独立生物学重复，也不自动合并。
- 6个HBB靶向样本：HBB WT、HBB IVS-110各3个重复。它们是小鼠细胞中的整合HBB位点实验，需要作者的靶向位点参考与相应计数方式，不直接套用全基因组SPI。
- ENA文件名含 `_subreads.fastq.gz`。仅据文件名/平台不能认证其为FLNC，也不能把FASTQ当成可重算CCS的原始SMRT movie或native BAM。本包保存归档序列，后续需核实实际CCS/read type、接头、方向、poly(A)过滤和SPI-PBLR入口契约。

来源：[GEO GSE144205](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144205)、[原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC8038867/)。原研究的常规比对使用mm10，HBB靶向数据使用自定义整合位点参考；未来选择统一参考版本时需另登记，当前未配置计算参考。

### scFLUENT-seq SPI

GSE278777中12个 **new RNA** 样本全部保留：10个小鼠ESC/Pol II AID相关测序样本，以及脾脏细胞和MEF各1个。该研究EU标记为10分钟；IAA处理和wash-off时长另列，不作为标记脉冲。

GEO全研究有16个样本；4个old RNA分支不在本次待下载新生RNA名单内。12个测序文库不是12个细胞，每个单细胞文库可包含多个细胞。

全部12个Run为双端，下载原始R1/R2并保留原始序列和header，不做bulk重命名、barcode/UMI剪切。scFLUENT计算前需确认真实barcode/UMI结构、read角色与白名单，按具体版本的专用Kit处理。

来源：[GEO GSE278777](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278777)。样本处理说明使用mm10/GENCODE vM22；这是研究参考线索，当前包没有自行创建该参考。

### scGRO-seq SR

GSE242176原表的82个候选中，76个样本标题和样本协议对应scGRO-seq；6个inAGTuC/AGTuC比较分支未纳入scGRO计算批次。它们没有被标为坏数据或从原表删除。

从76个scGRO候选选50个：优先保留32个PAIRED归档，再从SINGLE归档按实验批次覆盖选择18个。不是因深度/大小排除其他Run。

- 32个双端：原样保留R1/R2。
- 18个单端：原样保留ENA提供的单个FASTQ；单端元数据不能证明barcode被保留或被删除，仍需核对序列/原始header及plate/well对应关系。
- 96孔板barcodes与UMI由专用scGRO Kit处理；本包不把单端或双端的存在直接认作逐细胞SR可计算。

来源：[GEO GSE242176](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242176)、[原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC11222150/)。原研究使用mm10。完整所选Run与GSM、实验编号见`manifests/added_runs.tsv`。

## 新版下载器行为

- 原5类继续使用原来的prefetch → vdb-validate → fasterq-dump → FASTQ校验/文库输入准备流程。
- 新3类使用 **ENA原始归档FASTQ下载**，以官方文件大小、MD5、完整gzip读取、FASTQ格式和适用的双端ID/数量进行检查。保留原始文件字节和header，避免SRA通用导出时重新编号或忽略可能承载barcode的读段。
- 使用同步curl单文件断点下载；同时最多处理2个Run。curl禁用隐式配置、内部重试，外层重试最多3次，并等待当前传输结束。超长part、MD5不符和旧坏文件隔离保留；不使用旧大机的四连接分块代码。
- 失败不会终止其他Run；整批失败返回2。再次运行可复用指纹未变化且已有通过记录的文件。
- 新3类校验后生成 **`raw_fastq_manifest.json`** 与state，不生成虚假的FLNC/barcode准备完成标记；状态为`ARCHIVE_FASTQ_VALIDATED_KIT_PREPARATION_REQUIRED`。
- 本包不启动SPI/SR计算。3类仍需对应Kit的FLNC/单细胞adapter；Kit入口、参考和barcode参数未注册前，路由HELD。HBB靶向数据另有自定义参考限制。

## 输出位置

沿用用户指定目录：

```text
/srv/projects/CTS_Atlas/liuql/data/<Kit目录名>/<技术>/<RunID>/
  fastq/                 # 新3类的原始ENA FASTQ.gz，官方文件名保留
  sra/                   # 原5类的SRA缓存；新3类默认不额外下载SRA副本
  merged/                # 文库输入区；新3类等待专用adapter
  tmp/
  logs/ena_download.log
  stats/run_storage.json
  raw_fastq_manifest.json # 新3类原始文件验证凭据
  metadata.json
  routing.json
  state.json
```

`scFLUENT-seq SPI`和`scGRO-seq SR`对应安全目录名以`config/kit_registry.json`为准。包内`generated/`存逐Run结果、整批summary和文件大小统计。

## 运行命令

新服务器需Python3.9+、curl、可写数据目录。仅运行新增3类时不需要SRA Toolkit、Slurm、容器或参考文件；运行原5类/全部名单时可用包内安装器安装Toolkit。

将ZIP上传到 `/srv/projects/CTS_Atlas/liuql/`，解压到新的 `cts_download_merge_20261014` 目录。

### 仅下载这次新增的3类、76个Run

```bash
cd /srv/projects/CTS_Atlas/liuql/cts_download_merge_20261014
bash start.sh --check --only-added
mkdir -p logs
nohup bash start.sh --stage all --only-added --workers 2 \
  > logs/added_three_kits.out 2>&1 < /dev/null &
echo "PID=$!"
```

### 完整8类、339个Run

```bash
bash start.sh --install-tools --check
nohup bash start.sh --stage all --workers 2 \
  > logs/all_eight_kits.out 2>&1 < /dev/null &
echo "PID=$!"
```

原版已经下载的数据按原路径及已有校验状态复用。项目锁防止两个下载器同时写同一项目；若旧版仍运行，新版会拒绝启动，待旧版结束后再启动即可。不要为升级删掉已下载的数据。

查看进度：`tail -n 30 logs/added_three_kits.out`。结束后查看`generated/execution_summary.json`，逐Run信息在`generated/run_results.json`。大小可用`python3 tools/inventory_run_storage.py`刷新。

## 核对与测试

108个原始候选Run均核对了官方ENA记录，Run→SRX→GSM关系与GEO样本协议交叉匹配。原5类名单未变；新增76个Run没有与原名单重叠。证据/筛选摘要在`evidence/`，原工作簿未被修改。

在Ubuntu WSL上通过原SRA执行逻辑、启动入口、339个目录/278个样本路由、原始ENA文件下载/断点/校验/坏文件隔离的模拟测试。新3类未被标为最终可计算输入。测试不等于已在你的新服务器实测560.856GB下载。

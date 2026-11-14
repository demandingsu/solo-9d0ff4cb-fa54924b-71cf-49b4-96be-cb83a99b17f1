# 双端 FASTQ 质控 API

FastAPI 实现的成对 FASTQ 质量过滤服务：上传 `R1.fastq.gz` / `R2.fastq.gz`
与 Phred+33 质量阈值，后台流式处理，支持查询状态、停止任务、下载过滤后的
成对 gzip 结果。

## 快速启动

```bash
docker compose up --build
```

启动后：

* Swagger 文档：<http://localhost:8000/docs>
* OpenAPI JSON：<http://localhost:8000/openapi.json>
* 健康检查：<http://localhost:8000/health>

数据目录 `/data`（Compose 卷 `qc-data`）保存上传文件和结果；容器重启后
内存中的任务记录会失效，遗留目录会被自动清理。

不装 Docker 时也可本地运行：

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
QC_DATA_DIR=./data uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 接口一览

| 方法 | 路径 | 说明 |
|------|------|------|
| `POST` | `/jobs` | 上传 `r1`、`r2`（multipart `File`）和 `min_quality`，返回 `202` 与任务号 |
| `GET` | `/jobs/{job_id}` | 查询状态、进度、保留/丢弃对数、错误信息、下载链接 |
| `GET` | `/jobs` | 列出所有任务（新到旧） |
| `POST` | `/jobs/{job_id}/stop` | 请求停止；已完成/失败的任务返回 `409` |
| `GET` | `/jobs/{job_id}/download/r1` | 下载过滤后的 R1（仅成功任务） |
| `GET` | `/jobs/{job_id}/download/r2` | 下载过滤后的 R2（仅成功任务） |

`min_quality` 为 0..93 的整数（Phred+33 可表示的范围）。上传字段名固定为
`r1`、`r2`。

### 示例

```bash
# 提交任务
curl -X POST http://localhost:8000/jobs \
  -F 'r1=@samples/mixed_R1.fastq.gz' \
  -F 'r2=@samples/mixed_R2.fastq.gz' \
  -F 'min_quality=20'

# 查询
curl http://localhost:8000/jobs/<job_id>

# 停止
curl -X POST http://localhost:8000/jobs/<job_id>/stop

# 成功后可下载
curl -OJ http://localhost:8000/jobs/<job_id>/download/r1
curl -OJ http://localhost:8000/jobs/<job_id>/download/r2
```

### 状态字段

`status` 取值：`pending` → `running` → `succeeded` / `failed`，或
`running`/`pending` → `stopping` → `stopped`。

成功时：

```json
{
  "status": "succeeded",
  "total_pairs": 7,
  "kept_pairs": 4,
  "discarded_pairs": 3,
  "downloads": {
    "r1": "/jobs/<id>/download/r1",
    "r2": "/jobs/<id>/download/r2"
  },
  "error": null
}
```

失败时 `kept_pairs`/`discarded_pairs`/`downloads` 为 `null`，`error.code`
为下列之一：

| code | 含义 |
|------|------|
| `pair_name_mismatch` | 去掉 `/1`、`/2` 后缀后两端名称不同 |
| `invalid_name_suffix` | 名称未以 `/1`（R1）或 `/2`（R2）结尾 |
| `pair_count_mismatch` | 两端记录条数不同 |
| `seq_qual_length_mismatch` | 序列与质量行长度不同 |
| `malformed_record` | 记录不是完整四行、缺少 `@`/`+`、行超长等 |
| `invalid_quality_char` | 质量字符超出 Phred+33 可打印范围（ASCII 33..126） |
| `gzip_corrupt` | gzip 截断、CRC 错误、尾部垃圾、非 gzip 文件或空文件 |

## 处理规则

1. 每条读段恰好四行：`@名称`、序列、`+`、质量。
2. R1 名称必须以 `/1` 结尾，R2 必须以 `/2` 结尾；去掉后缀后逐对核对，
   必须完全一致。名称后的注释（空格分隔）会被忽略。
3. 序列与质量长度必须一致。
4. 质量按 Phred+33 解释；只要任一端出现低于阈值的碱基，整对丢弃。
5. 保留的读段按原顺序写出，输出名称沿用原始四行（仅补一个 `\n`）。
6. gzip 损坏、条数不符、名称不符、长度不符均导致任务 `failed`，不产生任何
   可下载文件。
7. 流式处理：逐对读取、逐对写出，不会整份载入内存；每 256 对更新一次进度。
   处理期间状态查询与停止请求由独立线程和事件立即响应。
8. 停止与完成互斥：一旦接受停止，任务不可能转为成功；完成后再次停止返回
   `409`。`stopping` 只在工作线程退出且临时文件（上传文件、`.part` 中间
   输出）清空后，才会变为 `stopped`。
9. 只有 `succeeded` 任务提供下载；失败或停止均返回 `409`，不会给出半成品。
   输出先在 `.part` 文件中完整写入，成功后才原子改名发布。

## 本地样例

`tools/make_samples.py` 可重新生成 `samples/` 下的样例：

| 文件 | 内容 |
|------|------|
| `good_R1/R2.fastq.gz` | 3 对全部通过的读段 |
| `bad_R1/R2.fastq.gz` | 3 对全部因低质量被丢弃 |
| `mixed_R1/R2.fastq.gz` | 7 对，前 3 与最后 1 对保留，中间 3 对丢弃（阈值 20） |
| `truncated_R1.fastq.gz` + `truncated_R2.fastq.gz` | R2 的 gzip 流末尾被截断，任务必定失败为 `gzip_corrupt` |

```bash
python tools/make_samples.py
```

## 测试

```bash
pip install -r requirements.txt pytest httpx
python -m pytest tests -q
```

测试覆盖过滤/顺序/阈值边界、各类强制失败、停止语义（含停止抢先于完成/
失败、清理时序）、处理期间查询响应性，以及流式读取行数上限。

## 说明与限制

* 单进程内存任务表，符合“单个 API 进程管理任务”的范围；不要用多 worker
  启动（`--workers 1`），否则任务查询会落到不同进程。Compose 默认单进程。
* 每个上传文件默认上限 20 GiB，可用 `QC_MAX_UPLOAD_BYTES` 调整。
* 后台线程为守护线程；进程关闭时会给所有任务发停止信号并等待清理。

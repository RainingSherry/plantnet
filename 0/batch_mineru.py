#!/usr/bin/env python3
"""
MinerU_Batch_Process.py
批量处理多个PDF分卷，每个分卷转换完成后合并Markdown
"""

import os
import sys
import subprocess
import time
from pathlib import Path

PDF_PARTS_DIR = "/home/luolie/biopipeline/dimension-reduction/plantnet/0/pdf_parts"
OUTPUT_BASE_DIR = "/home/luolie/biopipeline/dimension-reduction/plantnet/0/md_output"
MINERU_SCRIPT = "/home/luolie/.cursor/skills/mineru-pdf/scripts/mineru_upload.py"

def process_part(part_file, part_num):
    """处理单个PDF分卷"""
    print(f"\n{'='*60}")
    print(f"处理第 {part_num} 部分: {part_file.name}")
    print(f"{'='*60}\n")
    
    start_time = time.time()
    
    # 运行minerU上传脚本
    cmd = [
        "python3", MINERU_SCRIPT,
        str(part_file),
        "--language", "ch",
        "--model-version", "vlm",
        "--no-ocr-fallback"  # 跳过OCR加速处理
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=False)
        elapsed = time.time() - start_time
        print(f"\n第 {part_num} 部分完成，耗时: {elapsed/60:.1f} 分钟")
        return result.returncode == 0
    except Exception as e:
        print(f"处理失败: {e}")
        return False

def merge_markdown_files(output_dir, final_output):
    """合并所有Markdown文件"""
    print(f"\n{'='*60}")
    print("合并Markdown文件...")
    print(f"{'='*60}\n")
    
    md_files = sorted(Path(output_dir).glob("*.md"))
    print(f"找到 {len(md_files)} 个Markdown文件")
    
    merged_content = []
    merged_content.append("# 普林斯顿数学指南\n\n")
    merged_content.append(f"> 本文档由MinerU自动从PDF转换生成，共 {len(md_files)} 个分卷\n\n")
    merged_content.append("---\n\n")
    
    for i, md_file in enumerate(md_files, 1):
        print(f"读取: {md_file.name}")
        content = md_file.read_text(encoding='utf-8')
        merged_content.append(f"\n\n<!-- 第 {i} 部分: {md_file.stem} -->\n\n")
        merged_content.append(content)
    
    final_path = Path(final_output)
    final_path.write_text("".join(merged_content), encoding='utf-8')
    print(f"\n合并完成！输出文件: {final_path}")
    print(f"总大小: {final_path.stat().st_size / 1024 / 1024:.1f} MB")

def main():
    # 创建输出目录
    os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)
    
    # 获取所有PDF分卷
    pdf_files = sorted(Path(PDF_PARTS_DIR).glob("*.pdf"))
    print(f"找到 {len(pdf_files)} 个PDF分卷")
    
    if not pdf_files:
        print("未找到PDF文件！")
        sys.exit(1)
    
    # 依次处理每个分卷
    success_count = 0
    for i, pdf_file in enumerate(pdf_files, 1):
        if process_part(pdf_file, i):
            success_count += 1
        else:
            print(f"警告: 第 {i} 部分处理失败，继续下一部分...")
        
        # 每个分卷处理后休息一下，避免API限制
        if i < len(pdf_files):
            print("等待10秒后继续下一个分卷...")
            time.sleep(10)
    
    print(f"\n处理完成！成功: {success_count}/{len(pdf_files)}")
    
    # 合并Markdown
    if success_count > 0:
        final_output = Path(OUTPUT_BASE_DIR) / "普林斯顿数学指南_合并版.md"
        merge_markdown_files(OUTPUT_BASE_DIR, final_output)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Princeton_Mathematical_Guide_PDF_Splitter.py
将大型PDF分割成多个小批次，每个批次约150页
"""

import fitz
import os
from pathlib import Path

def split_pdf(input_pdf, output_dir, pages_per_batch=150):
    """分割PDF为多个批次"""
    doc = fitz.open(input_pdf)
    total_pages = len(doc)
    num_batches = (total_pages + pages_per_batch - 1) // pages_per_batch
    
    print(f"总页数: {total_pages}")
    print(f"每批页数: {pages_per_batch}")
    print(f"总共 {num_batches} 批")
    
    os.makedirs(output_dir, exist_ok=True)
    
    batch_files = []
    for i in range(num_batches):
        start_page = i * pages_per_batch + 1  # 1-indexed
        end_page = min((i + 1) * pages_per_batch, total_pages)
        
        output_pdf = os.path.join(output_dir, f"part_{i+1:03d}_pages_{start_page}-{end_page}.pdf")
        
        # 创建新文档并添加页面
        new_doc = fitz.open()
        for page_num in range(start_page - 1, end_page):  # 0-indexed
            new_doc.insert_pdf(doc, from_page=page_num, to_page=page_num)
        
        new_doc.save(output_pdf)
        new_doc.close()
        
        batch_files.append(output_pdf)
        print(f"已创建: {output_pdf} (页 {start_page}-{end_page})")
    
    doc.close()
    return batch_files

if __name__ == "__main__":
    input_pdf = "/home/luolie/biopipeline/dimension-reduction/plantnet/0/普林斯顿数学指南.pdf"
    output_dir = "/home/luolie/biopipeline/dimension-reduction/plantnet/0/pdf_parts"
    
    batch_files = split_pdf(input_pdf, output_dir, pages_per_batch=150)
    print(f"\n分割完成！共 {len(batch_files)} 个文件")

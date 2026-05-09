output_file="trj_reco_0921.xyz"  # 合并后的输出文件名

# 循环遍历目录A1到A10，将所有txt文件路径添加到一个数组中
    file_list=()  
    for i in $(seq 1 10 ); do
        file_list+=("A$i/trj_reco_0921.xyz")
    done

    # 清空输出文件
    > "$output_file"

    # 循环遍历文件列表，逐个将内容追加到输出文件
    for file in "${file_list[@]}"
    do
        cat "$file" >> "$output_file"
    done

    echo "文件合并完成！结果保存在 $output_file 文件中。"


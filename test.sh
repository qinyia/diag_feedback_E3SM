aa=(a b c)
echo ${aa[@]}
bb=(dd ee f)
echo ${bb[@]}
cc=( "${aa[@]}" "${bb[@]}")
echo ${cc[@]}

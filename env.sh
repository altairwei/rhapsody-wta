export PATH="$PWD/bin:$PWD/scripts:$PATH"

function proj_dir {
    mkdir -p data/annotations
    mkdir -p data/reference_indexes
    mkdir -p data/reference_sequences

    mkdir -p logs
    mkdir -p tmp

    mkdir -p results
}

function bytes_readable() {
	local size=${1:-"$(cat)"} factor="KMGTEPZY" scale="scale=2"
	if (( ${size} < 1024 )); then
		echo "${size} bytes"
		return 0
	else
		size=$(echo "${scale}; ${size}/1024" | bc)
	fi
	while (( $(echo "${size} >= 1024" | bc -l) && ${#factor} > 1 )); do
		size=$(echo "${scale}; ${size}/1024" | bc)
		factor=${factor:1}
	done
	echo "${size} ${factor:0:1}iB"
}

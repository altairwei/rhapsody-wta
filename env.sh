die()
{
        local _ret=$2
        test -n "$_ret" || _ret=1
        test "$_PRINT_HELP" = yes && print_help >&2
        echo "$1" >&2
        exit ${_ret}
}


SCRIPT_DIR="$(cd "$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")" && pwd)" || die "Couldn't determine the script's running directory, which probably matters, bailing out" 2

export PATH="${SCRIPT_DIR}/bin:${SCRIPT_DIR}/scripts:$PATH"
export PYTHONPATH="${SCRIPT_DIR}/docker/rhapsody_1.8/mist/src:$PYTHONPATH"

function proj_dir {
	mkdir -p ${SCRIPT_DIR}/data/reads
    mkdir -p ${SCRIPT_DIR}/data/annotations
    mkdir -p ${SCRIPT_DIR}/data/reference_indexes
    mkdir -p ${SCRIPT_DIR}/data/reference_sequences

    mkdir -p ${SCRIPT_DIR}/logs
    mkdir -p ${SCRIPT_DIR}/tmp

    mkdir -p ${SCRIPT_DIR}/results
}

function bytes_readable {
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
	echo "${size} ${factor:0:1}B"
}

function clean_tmp {
	rm -r ${SCRIPT_DIR}/tmp*
	rm -r ${SCRIPT_DIR}/node-*
	rm -r ${SCRIPT_DIR}/*-cleanup-arena-members

	mkdir -p ${SCRIPT_DIR}/tmp/cwltool
	mkdir -p ${SCRIPT_DIR}/tmp/toil
}
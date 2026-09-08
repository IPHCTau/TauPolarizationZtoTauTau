#!/usr/bin/env bash

# Usage:
#   source slurm_setup.sh
# Run this script from a fresh shell in the analysis repository.

# Use a locale available on IPHC login and worker nodes.
unset LANGUAGE LC_ALL LC_CTYPE
export LANG=C.utf8
export LC_ALL=C.utf8
export LC_CTYPE=C.utf8

# Load Columnflow and all paths configured in .setups/zttenv.sh.
source setup.sh zttenv || return 1

# Restore the worker-compatible locale after environment activation.
unset LANGUAGE
export LANG=C.utf8
export LC_ALL=C.utf8
export LC_CTYPE=C.utf8

# The same filesystem has different physical mount points:
#
#   login physical path : /mnt/slurm/shared
#   common alias        : /slurm/shared
#   worker physical path: /grid-mnt/shared
#
# Only /slurm/shared works with the same spelling on login and workers.
canonicalize_value() {
    local current_value="${1:-}"
    current_value="${current_value//\/mnt\/slurm\/shared/\/slurm\/shared}"
    printf '%s' "$current_value"
}

canonicalize_environment_paths() {
    local variable_name
    local variable_value

    # Canonicalize every exported Columnflow variable.
    while IFS='=' read -r variable_name variable_value; do
        [[ "$variable_name" == CF_* ]] || continue
        variable_value="$(canonicalize_value "$variable_value")"
        export "$variable_name=$variable_value"
    done < <(env)

    # Canonicalize environment variables that can contain conda/venv paths.
    local path_variables=(
        MAMBA_ROOT_PREFIX
        CONDA_PREFIX
        VIRTUAL_ENV
        PATH
        PYTHONPATH
        LD_LIBRARY_PATH
        LIBRARY_PATH
        CPATH
        CMAKE_PREFIX_PATH
        PKG_CONFIG_PATH
    )

    for variable_name in "${path_variables[@]}"; do
        if [[ -n "${!variable_name:-}" ]]; then
            variable_value="$(canonicalize_value "${!variable_name}")"
            export "$variable_name=$variable_value"
        fi
    done
}

canonicalize_environment_paths

# Validate variables needed by local setup and remote Slurm jobs.
required_variables=(
    CF_BASE
    CF_REPO_BASE
    CF_DATA
    CF_SOFTWARE_BASE
    CF_CONDA_BASE
    CF_VENV_BASE
    CF_JOB_BASE
    CF_STORE_LOCAL
)

for variable_name in "${required_variables[@]}"; do
    if [[ -z "${!variable_name:-}" ]]; then
        echo "ERROR: $variable_name is unset or empty"
        return 1
    fi
done

# Repair paths generated inside micromamba.sh on the login node.
fix_micromamba_worker_path() {
    local mamba_init_file="$CF_CONDA_BASE/etc/profile.d/micromamba.sh"

    if [[ ! -f "$mamba_init_file" ]]; then
        echo "ERROR: micromamba initialization file does not exist:"
        echo "  $mamba_init_file"
        return 1
    fi

    # Fix fully expanded and dynamically constructed login-only paths.
    sed -i \
        -e 's|/mnt/slurm/shared|/slurm/shared|g' \
        -e 's|/mnt/grid${MAMBA_ROOT_PREFIX}|${MAMBA_ROOT_PREFIX}|g' \
        -e 's|/mnt${MAMBA_ROOT_PREFIX}|${MAMBA_ROOT_PREFIX}|g' \
        "$mamba_init_file" || return 1

    if grep -qE \
        '/mnt/slurm/shared|/mnt(/grid)?\$\{MAMBA_ROOT_PREFIX\}' \
        "$mamba_init_file"
    then
        echo "ERROR: micromamba.sh contains worker-incompatible paths:"
        grep -nE \
            '/mnt/slurm/shared|/mnt(/grid)?\$\{MAMBA_ROOT_PREFIX\}' \
            "$mamba_init_file"
        return 1
    fi

    export MAMBA_ROOT_PREFIX="$CF_CONDA_BASE"

    echo "micromamba.sh uses the worker-safe path:"
    sed -n '5,8p' "$mamba_init_file"
}

fix_micromamba_worker_path || return 1

# Validate the shared software installation.
required_software_paths=(
    "$CF_SOFTWARE_BASE"
    "$CF_CONDA_BASE"
    "$CF_CONDA_BASE/bin/micromamba"
    "$CF_VENV_BASE"
)

for required_path in "${required_software_paths[@]}"; do
    if [[ ! -e "$required_path" ]]; then
        echo "ERROR: required software path does not exist:"
        echo "  $required_path"
        return 1
    fi
done

if [[ ! -x "$CF_CONDA_BASE/bin/micromamba" ]]; then
    echo "ERROR: micromamba is not executable:"
    echo "  $CF_CONDA_BASE/bin/micromamba"
    return 1
fi

# Create and validate job, output and fallback-status directories.
for directory in "$CF_JOB_BASE" "$CF_STORE_LOCAL"; do
    if [[ ! -d "$directory" ]]; then
        mkdir -p "$directory" || {
            echo "ERROR: could not create directory:"
            echo "  $directory"
            return 1
        }
        echo "Created directory: $directory"
    fi

    if [[ ! -r "$directory" || ! -w "$directory" ]]; then
        echo "ERROR: directory is not readable and writable:"
        echo "  $directory"
        return 1
    fi
done

export CF_SLURM_STATUS_DIR="${CF_JOB_BASE}/slurm_status"

if [[ ! -d "$CF_SLURM_STATUS_DIR" ]]; then
    mkdir -p "$CF_SLURM_STATUS_DIR" || {
        echo "ERROR: could not create Slurm status directory:"
        echo "  $CF_SLURM_STATUS_DIR"
        return 1
    }
    echo "Created Slurm status directory: $CF_SLURM_STATUS_DIR"
else
    echo "Slurm status directory exists: $CF_SLURM_STATUS_DIR"
fi

if [[ ! -r "$CF_SLURM_STATUS_DIR" || ! -w "$CF_SLURM_STATUS_DIR" ]]; then
    echo "ERROR: Slurm status directory is not readable and writable:"
    echo "  $CF_SLURM_STATUS_DIR"
    return 1
fi

# Required by tasks that access files through WLCG protocols.
export CF_WLCG_TOOLS="$CF_BASE/modules/law/src/law/contrib/wlcg/scripts/law_wlcg_tools.sh"

if [[ ! -r "$CF_WLCG_TOOLS" ]]; then
    echo "ERROR: CF_WLCG_TOOLS is not readable:"
    echo "  $CF_WLCG_TOOLS"
    return 1
fi

# Perform one final repair and validation immediately before submission.
canonicalize_environment_paths
fix_micromamba_worker_path || return 1

mamba_init="$CF_CONDA_BASE/etc/profile.d/micromamba.sh"

if grep -qE \
    '/mnt/slurm/shared|/mnt(/grid)?\$\{MAMBA_ROOT_PREFIX\}' \
    "$mamba_init"
then
    echo "ERROR: worker-incompatible micromamba path remains"
    return 1
fi

for variable_name in \
    CF_DATA \
    CF_SOFTWARE_BASE \
    CF_CONDA_BASE \
    CF_VENV_BASE \
    CF_JOB_BASE \
    CF_STORE_LOCAL \
    MAMBA_ROOT_PREFIX
do
    if [[ "${!variable_name}" == *"/mnt/slurm/shared"* ]]; then
        echo "ERROR: $variable_name contains /mnt/slurm/shared:"
        echo "  ${!variable_name}"
        return 1
    fi
done

echo
echo "Slurm environment ready"
echo "CF_REPO_BASE     : $CF_REPO_BASE"
echo "CF_DATA          : $CF_DATA"
echo "CF_SOFTWARE_BASE : $CF_SOFTWARE_BASE"
echo "CF_CONDA_BASE    : $CF_CONDA_BASE"
echo "CF_VENV_BASE     : $CF_VENV_BASE"
echo "CF_JOB_BASE      : $CF_JOB_BASE"
echo "CF_STORE_LOCAL   : $CF_STORE_LOCAL"
echo "Slurm status dir : $CF_SLURM_STATUS_DIR"
echo "MAMBA_ROOT_PREFIX: $MAMBA_ROOT_PREFIX"
echo "Sandbox handling : delegated to Slurm task dependencies"

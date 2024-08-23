import subprocess
import sys
import os

def install_from_local_git(repo_path, package_name):
    # Check if the repository exists locally
    if not os.path.exists(repo_path):
        print(f"Cloning the repository to {repo_path}...")
        git_clone_command = ["git", "clone", repo_path]
        subprocess.check_call(git_clone_command)
    else:
        print(f"Repository found at {repo_path}, pulling latest changes...")
        git_pull_command = ["git", "-C", repo_path, "pull"]
        subprocess.check_call(git_pull_command)

    # Install the package using pip from the local repository
    print(f"Installing {package_name} from the local repository...")
    pip_install_command = [sys.executable, "-m", "pip", "install", repo_path]
    subprocess.check_call(pip_install_command)

def import_or_install(package_name, repo_path):
    try:
        __import__(package_name)
        print(f"{package_name} is already installed.")
    except ImportError:
        print(f"{package_name} not found, installing...")
        install_from_local_git(repo_path, package_name)

# Example usage:
# Define the path to the local Git repository
local_repo_path = "./freesasa-repo"

# Attempt to import or install the freesasa package
import_or_install("freesasa", local_repo_path)

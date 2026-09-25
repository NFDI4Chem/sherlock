#!/bin/bash

# Define variables
PROJECT_DIR="/mnt/data/sherlock"
COMPOSE_FILE="docker-compose.yml"
SHERLOCK_IMAGE="nfdi4chem/sherlock:latest"
LOG_FILE="/var/log/sherlock-deploy.log"
LOG_OWNER="${SUDO_USER:-$(whoami)}"

# Create log file if it doesn't exist
if [ ! -f "$LOG_FILE" ]; then
    sudo touch "$LOG_FILE"
    sudo chmod 644 "$LOG_FILE"
    sudo chown "$LOG_OWNER":"$LOG_OWNER" "$LOG_FILE"
fi

# Unified logging function
log_message() {
    echo "$1"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$LOG_FILE"
}

# === Start of script ===
log_message "=========================================="
log_message "Starting Sherlock Deployment Script"
log_message "=========================================="

# Change to project directory to ensure paths resolve correctly
cd "$PROJECT_DIR" || {
    log_message "Failed to change to directory $PROJECT_DIR"
    exit 1
}
log_message "Working directory: $(pwd)"

# === Functions ===

# Cleanup function
cleanup() {
    log_message "Cleaning up dangling images..."
    docker image prune -f >/dev/null 2>&1 || true
    log_message "Cleanup completed"
}

# Deploy a service by pulling latest image and recreating container if updated
deploy_service() {
    local service_name=$1
    local image=$2

    log_message "Starting deployment for service: $service_name"
    log_message "Checking for new image: $image"

    # Pull the latest image
    if [ "$(docker pull "$image" | grep -c "Status: Image is up to date")" -eq 0 ]; then
        log_message "New image detected for $service_name"
        log_message "Recreating container with updated image..."
        docker compose -f "$COMPOSE_FILE" up -d --force-recreate --no-deps "$service_name"
        log_message "Deployment of $service_name completed successfully"
    else
        log_message "Image for $service_name is up to date. Skipping deployment."
    fi
}

# Main deployment process
main() {
    log_message "----------------------------------------"
    log_message "Deploying Sherlock Service"
    log_message "----------------------------------------"
    deploy_service "sherlock" "$SHERLOCK_IMAGE"

    log_message ""
    cleanup

    log_message ""
    log_message "=========================================="
    log_message "All Deployments Completed Successfully!"
    log_message "=========================================="
}

# Execute main deployment
main

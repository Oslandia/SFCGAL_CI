/**
 * SFCGAL WebAssembly Module Loader
 * Shared utility to load and initialize SFCGAL in all examples
 */

let sfcgalInstance = null;

export async function loadSFCGAL() {
    if (sfcgalInstance) {
        return sfcgalInstance;
    }

    try {
        const module = await import('./sfcgal.js');
        const Module = await module.default();

        if (!Module.SFCGAL) {
            throw new Error('SFCGAL class not found in module');
        }

        sfcgalInstance = new Module.SFCGAL();
        sfcgalInstance.initialize();

        return sfcgalInstance;
    } catch (error) {
        console.error('Failed to load SFCGAL:', error);
        throw error;
    }
}

export function showStatus(elementId, message, type = 'info') {
    const statusEl = document.getElementById(elementId);
    if (!statusEl) return;

    const classes = {
        'info': 'status-info',
        'success': 'status-success',
        'error': 'status-error',
        'warning': 'status-warning'
    };

    statusEl.className = `status ${classes[type] || classes.info}`;
    statusEl.textContent = message;
    statusEl.style.display = 'block';
}

export function hideStatus(elementId, delay = 0) {
    const statusEl = document.getElementById(elementId);
    if (!statusEl) return;

    if (delay > 0) {
        setTimeout(() => {
            statusEl.style.display = 'none';
        }, delay);
    } else {
        statusEl.style.display = 'none';
    }
}

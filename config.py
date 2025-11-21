"""
Centralized configuration for the research-agent system.
All configurable values should be defined here.
"""
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# =============================================================================
# LLM CONFIGURATION
# =============================================================================

# Provider settings
LLM_DEFAULT_PROVIDER = os.getenv('LLM_DEFAULT_PROVIDER', 'mock')

# Model IDs
CLAUDE_MODEL = os.getenv('CLAUDE_MODEL', 'claude-sonnet-4-5-20250514')
GPT4_MODEL = os.getenv('GPT4_MODEL', 'gpt-4-turbo')

# Generation defaults
DEFAULT_TEMPERATURE = float(os.getenv('DEFAULT_TEMPERATURE', '0.7'))
DEFAULT_MAX_TOKENS = int(os.getenv('DEFAULT_MAX_TOKENS', '4096'))

# Retry settings
LLM_MAX_RETRIES = int(os.getenv('LLM_MAX_RETRIES', '3'))
LLM_RETRY_DELAY = float(os.getenv('LLM_RETRY_DELAY', '1.0'))

# =============================================================================
# STORAGE CONFIGURATION
# =============================================================================

# Directories
PROJECT_DATA_DIR = os.getenv('PROJECT_DATA_DIR', 'data/projects')
OUTPUTS_DIR = os.getenv('OUTPUTS_DIR', 'outputs')
LOG_DIR = os.getenv('LOG_DIR', 'logs')

# =============================================================================
# DASHBOARD CONFIGURATION
# =============================================================================

# UI settings
QUICK_AGENT_LIMIT = int(os.getenv('QUICK_AGENT_LIMIT', '5'))
PROJECT_LOAD_TIMEOUT = int(os.getenv('PROJECT_LOAD_TIMEOUT', '30000'))  # ms

# =============================================================================
# COST ESTIMATION (per 1000 tokens)
# =============================================================================

ANTHROPIC_INPUT_COST = float(os.getenv('ANTHROPIC_INPUT_COST', '0.003'))
ANTHROPIC_OUTPUT_COST = float(os.getenv('ANTHROPIC_OUTPUT_COST', '0.015'))
GPT4_INPUT_COST = float(os.getenv('GPT4_INPUT_COST', '0.01'))
GPT4_OUTPUT_COST = float(os.getenv('GPT4_OUTPUT_COST', '0.03'))

# =============================================================================
# API KEY VALIDATION
# =============================================================================

ANTHROPIC_KEY_PREFIX = 'sk-ant-'
OPENAI_KEY_PREFIX = 'sk-'

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

LOG_LEVEL = os.getenv('LOG_LEVEL', 'INFO')
LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

# =============================================================================
# AGENT CONFIGURATION
# =============================================================================

# Temperature overrides for specific agents
AGENT_TEMPERATURES = {
    'orchestrator': 0.5,      # More deterministic planning
    'research_prd': 0.7,      # Balanced creativity
    'literature': 0.6,        # Slightly more focused
    'framework': 0.7,
    'data_qa': 0.3,           # Very deterministic
    'eda': 0.5,
    'modeling': 0.4,          # More deterministic
    'figures': 0.6,
    'writer': 0.8,            # More creative
    'references': 0.2,        # Very deterministic
    'reviewer': 0.5,
    'submission': 0.4,
}

# Max tokens for specific agents
AGENT_MAX_TOKENS = {
    'orchestrator': 3000,
    'research_prd': 8000,
    'literature': 6000,
    'framework': 4000,
    'data_qa': 4000,
    'eda': 5000,
    'modeling': 5000,
    'figures': 4000,
    'writer': 8000,
    'references': 3000,
    'reviewer': 5000,
    'submission': 4000,
}


def get_agent_temperature(agent_type: str) -> float:
    """Get temperature for a specific agent type."""
    return AGENT_TEMPERATURES.get(agent_type, DEFAULT_TEMPERATURE)


def get_agent_max_tokens(agent_type: str) -> int:
    """Get max tokens for a specific agent type."""
    return AGENT_MAX_TOKENS.get(agent_type, DEFAULT_MAX_TOKENS)


def validate_api_key(provider: str, key: str) -> tuple[bool, str]:
    """Validate API key format for a provider."""
    if not key or not key.strip():
        return False, "API key cannot be empty"

    if provider == 'anthropic':
        if not key.startswith(ANTHROPIC_KEY_PREFIX):
            return False, f"Anthropic API key must start with '{ANTHROPIC_KEY_PREFIX}'"
    elif provider == 'openai':
        if not key.startswith(OPENAI_KEY_PREFIX):
            return False, f"OpenAI API key must start with '{OPENAI_KEY_PREFIX}'"

    return True, ""

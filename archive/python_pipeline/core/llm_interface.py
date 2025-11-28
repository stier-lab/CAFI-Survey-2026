"""
LLM Interface Module - Handles integration with various LLM providers
"""
from typing import Dict, List, Any, Optional
from abc import ABC, abstractmethod
import os
import time
import logging

# Set up logging
logger = logging.getLogger(__name__)

# Try to import config, fall back to defaults
try:
    from config import (
        CLAUDE_MODEL, GPT4_MODEL, DEFAULT_TEMPERATURE, DEFAULT_MAX_TOKENS,
        LLM_MAX_RETRIES, LLM_RETRY_DELAY
    )
except ImportError:
    CLAUDE_MODEL = 'claude-sonnet-4-20250514'
    GPT4_MODEL = 'gpt-4-turbo'
    DEFAULT_TEMPERATURE = 0.7
    DEFAULT_MAX_TOKENS = 4096
    LLM_MAX_RETRIES = 3
    LLM_RETRY_DELAY = 1.0


class LLMProvider(ABC):
    """Base class for LLM providers."""

    @abstractmethod
    def generate(self, system_prompt: str, user_prompt: str,
                max_tokens: int = 4096, temperature: float = 0.7) -> str:
        """Generate text using the LLM."""
        pass


class AnthropicProvider(LLMProvider):
    """Anthropic Claude provider with retry logic and error handling."""

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")

        if not self.api_key:
            raise ValueError("Anthropic API key not provided. Set ANTHROPIC_API_KEY environment variable.")

        try:
            import anthropic
            self.client = anthropic.Anthropic(api_key=self.api_key)
            self.anthropic = anthropic  # Store module for exception handling
        except ImportError:
            raise ImportError("anthropic package not installed. Run: pip install anthropic")

        # Usage tracking
        self.usage_stats = {
            'total_input_tokens': 0,
            'total_output_tokens': 0,
            'total_calls': 0,
            'total_errors': 0
        }

    def generate(self, system_prompt: str, user_prompt: str,
                max_tokens: int = 4096, temperature: float = 0.7) -> str:
        """Generate text using Claude with retry logic."""

        for attempt in range(LLM_MAX_RETRIES):
            try:
                message = self.client.messages.create(
                    model=CLAUDE_MODEL,
                    max_tokens=min(max_tokens, 8192),  # Cap at API limit
                    temperature=max(0.0, min(1.0, temperature)),  # Validate range
                    system=system_prompt,
                    messages=[{"role": "user", "content": user_prompt}]
                )

                # Track usage
                if hasattr(message, 'usage'):
                    self.usage_stats['total_input_tokens'] += message.usage.input_tokens
                    self.usage_stats['total_output_tokens'] += message.usage.output_tokens
                self.usage_stats['total_calls'] += 1

                if not message.content:
                    raise ValueError("Empty response from Claude API")

                return message.content[0].text

            except self.anthropic.RateLimitError as e:
                self.usage_stats['total_errors'] += 1
                if attempt < LLM_MAX_RETRIES - 1:
                    wait_time = LLM_RETRY_DELAY * (2 ** attempt)
                    logger.warning(f"Rate limited, waiting {wait_time}s before retry...")
                    time.sleep(wait_time)
                    continue
                logger.error(f"Rate limit exceeded after {LLM_MAX_RETRIES} retries")
                raise

            except self.anthropic.APIError as e:
                self.usage_stats['total_errors'] += 1
                if attempt < LLM_MAX_RETRIES - 1 and hasattr(e, 'status_code') and e.status_code >= 500:
                    logger.warning(f"API error {e.status_code}, retrying...")
                    time.sleep(LLM_RETRY_DELAY)
                    continue
                logger.error(f"Anthropic API error: {e}")
                raise

            except Exception as e:
                self.usage_stats['total_errors'] += 1
                logger.error(f"Unexpected error in generate: {e}")
                raise

    def get_usage_stats(self) -> Dict[str, int]:
        """Get token usage statistics."""
        return self.usage_stats.copy()


class OpenAIProvider(LLMProvider):
    """OpenAI GPT provider."""

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key or os.environ.get("OPENAI_API_KEY")

        if not self.api_key:
            raise ValueError("OpenAI API key not provided")

        try:
            from openai import OpenAI
            self.client = OpenAI(api_key=self.api_key)
        except ImportError:
            raise ImportError("openai package not installed. Run: pip install openai")

    def generate(self, system_prompt: str, user_prompt: str,
                max_tokens: int = 4096, temperature: float = 0.7) -> str:
        """Generate text using GPT."""
        response = self.client.chat.completions.create(
            model="gpt-4",
            max_tokens=max_tokens,
            temperature=temperature,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt}
            ]
        )

        return response.choices[0].message.content


class MockProvider(LLMProvider):
    """Mock provider for testing without API calls."""

    def generate(self, system_prompt: str, user_prompt: str,
                max_tokens: int = 4096, temperature: float = 0.7) -> str:
        """Return template response for testing."""
        return f"[Mock LLM Response]\n\nSystem: {system_prompt[:50]}...\n\nUser: {user_prompt[:50]}...\n\n[This is a placeholder response. To get real AI-generated content, configure an LLM provider.]"


class LLMInterface:
    """Main interface for LLM interactions."""

    def __init__(self, provider: str = "mock", api_key: Optional[str] = None):
        """
        Initialize LLM interface.

        Args:
            provider: "anthropic", "openai", or "mock"
            api_key: API key for the provider (or set via environment variable)
        """
        self.provider_name = provider

        if provider == "anthropic":
            self.provider = AnthropicProvider(api_key)
        elif provider == "openai":
            self.provider = OpenAIProvider(api_key)
        elif provider == "mock":
            self.provider = MockProvider()
        else:
            raise ValueError(f"Unknown provider: {provider}")

    def generate(self, system_prompt: str, user_prompt: str,
                max_tokens: int = 4096, temperature: float = 0.7) -> str:
        """
        Generate text using the configured LLM provider.

        Args:
            system_prompt: System instructions for the LLM
            user_prompt: User prompt/query
            max_tokens: Maximum tokens to generate
            temperature: Temperature for generation (0.0-1.0)

        Returns:
            Generated text
        """
        return self.provider.generate(system_prompt, user_prompt, max_tokens, temperature)

    def batch_generate(self, prompts: List[Dict[str, str]],
                      max_tokens: int = 4096, temperature: float = 0.7) -> List[str]:
        """
        Generate multiple responses in batch.

        Args:
            prompts: List of dicts with 'system_prompt' and 'user_prompt'
            max_tokens: Maximum tokens per response
            temperature: Temperature for generation

        Returns:
            List of generated texts
        """
        results = []
        for prompt in prompts:
            result = self.generate(
                prompt['system_prompt'],
                prompt['user_prompt'],
                max_tokens,
                temperature
            )
            results.append(result)

        return results


# Global LLM instance (can be configured once for all agents)
_global_llm: Optional[LLMInterface] = None


def configure_llm(provider: str = "mock", api_key: Optional[str] = None):
    """Configure the global LLM instance."""
    global _global_llm
    _global_llm = LLMInterface(provider, api_key)


def get_llm() -> LLMInterface:
    """Get the global LLM instance."""
    global _global_llm
    if _global_llm is None:
        # Default to mock provider if not configured
        _global_llm = LLMInterface("mock")
    return _global_llm


# Auto-detect and configure based on environment variables
def auto_configure_llm():
    """Automatically configure LLM based on available API keys."""
    if os.environ.get("ANTHROPIC_API_KEY"):
        try:
            configure_llm("anthropic")
            return "anthropic"
        except Exception:
            pass

    if os.environ.get("OPENAI_API_KEY"):
        try:
            configure_llm("openai")
            return "openai"
        except Exception:
            pass

    # Fall back to mock
    configure_llm("mock")
    return "mock"


# Auto-configure on module import
_configured_provider = auto_configure_llm()

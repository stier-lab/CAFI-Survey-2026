"""
Base Agent Class for Research Workflow System
"""
from abc import ABC, abstractmethod
from typing import Dict, List, Any, Optional
from datetime import datetime
import json
import os
import sys

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class BaseAgent(ABC):
    """Base class for all research workflow agents."""

    def __init__(self, name: str, role: str, mission: str):
        self.name = name
        self.role = role
        self.mission = mission
        self.execution_history: List[Dict[str, Any]] = []

    @abstractmethod
    def get_system_prompt(self) -> str:
        """Return the system prompt for this agent."""
        pass

    @abstractmethod
    def get_inputs(self) -> List[str]:
        """Return list of required input documents/data."""
        pass

    @abstractmethod
    def get_outputs(self) -> List[str]:
        """Return list of expected output artifacts."""
        pass

    @abstractmethod
    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute the agent's core functionality.

        Args:
            inputs: Dictionary of input documents/data
            context: Additional context including project brief, previous outputs, etc.

        Returns:
            Dictionary containing outputs and metadata
        """
        pass

    def validate_inputs(self, inputs: Dict[str, Any]) -> tuple[bool, str]:
        """
        Validate that all required inputs are provided and non-empty.

        Returns:
            Tuple of (is_valid, error_message)
        """
        required = self.get_inputs()
        missing = []
        empty = []

        for req in required:
            if req not in inputs:
                missing.append(req)
            elif inputs[req] is None or (isinstance(inputs[req], str) and not inputs[req].strip()):
                empty.append(req)
            elif isinstance(inputs[req], str) and inputs[req].startswith('[') and inputs[req].endswith('not provided]'):
                empty.append(req)

        if missing:
            return False, f"Missing required inputs: {', '.join(missing)}"
        if empty:
            return False, f"Empty values for required inputs: {', '.join(empty)}"

        return True, ""

    def validate_inputs_legacy(self, inputs: Dict[str, Any]) -> bool:
        """Legacy validation method for backwards compatibility."""
        is_valid, _ = self.validate_inputs(inputs)
        return is_valid

    def log_execution(self, inputs: Dict[str, Any], outputs: Dict[str, Any],
                     status: str, error: Optional[str] = None):
        """Log execution details for tracking."""
        self.execution_history.append({
            'timestamp': datetime.now().isoformat(),
            'inputs': list(inputs.keys()),
            'outputs': list(outputs.keys()) if outputs else [],
            'status': status,
            'error': error
        })

    def get_metadata(self) -> Dict[str, Any]:
        """Return agent metadata."""
        return {
            'name': self.name,
            'role': self.role,
            'mission': self.mission,
            'required_inputs': self.get_inputs(),
            'expected_outputs': self.get_outputs(),
            'execution_count': len(self.execution_history)
        }

    def format_prompt(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> str:
        """Format the complete prompt for the LLM including system prompt and inputs."""
        prompt = "=== PROJECT CONTEXT ===\n\n"

        if 'project_brief' in context:
            prompt += f"PROJECT BRIEF:\n{context['project_brief']}\n\n"

        prompt += "=== INPUTS ===\n\n"
        for key, value in inputs.items():
            prompt += f"{key.upper()}:\n{value}\n\n"

        prompt += "\n=== TASK ===\n\n"
        prompt += f"Generate the outputs as specified in the system prompt for the {self.name}.\n"
        prompt += "Provide comprehensive, well-structured, and academically rigorous content.\n"

        return prompt

    def call_llm(self, inputs: Dict[str, Any], context: Dict[str, Any],
                temperature: float = 0.7, max_tokens: int = 4096) -> str:
        """
        Call LLM with formatted prompts.

        Args:
            inputs: Input dictionary
            context: Context dictionary
            temperature: LLM temperature
            max_tokens: Maximum tokens to generate

        Returns:
            LLM generated text
        """
        try:
            from core.llm_interface import get_llm

            llm = get_llm()
            system_prompt = self.get_system_prompt()
            user_prompt = self.format_prompt(inputs, context)

            response = llm.generate(system_prompt, user_prompt, max_tokens, temperature)
            return response

        except ImportError:
            # Fallback if LLM interface not available
            return "[LLM not configured - using template response]"
        except Exception as e:
            raise RuntimeError(f"LLM call failed: {str(e)}")

"""
Core modules for Research Workflow System
"""
from .base_agent import BaseAgent
from .project_state import ProjectState, ProjectStage, TaskStatus

__all__ = [
    'BaseAgent',
    'ProjectState',
    'ProjectStage',
    'TaskStatus',
]

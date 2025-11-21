"""
Research Workflow Agents Package for CAFI Survey Analysis
"""
from .orchestrator_agent import OrchestratorAgent
from .research_prd_agent import ResearchPRDAgent
from .all_agents import (
    DeepLiteratureAgent,
    ConceptualFrameworkAgent,
    DataQAAgent,
    EDAAgent,
    ModelingAgent,
    FigureFactoryAgent,
    ScientificWriterAgent,
    ReferenceAgent,
    ReviewerAgent,
    SubmissionAgent,
    get_agent,
    list_agents,
    AGENT_REGISTRY
)

__all__ = [
    'OrchestratorAgent',
    'ResearchPRDAgent',
    'DeepLiteratureAgent',
    'ConceptualFrameworkAgent',
    'DataQAAgent',
    'EDAAgent',
    'ModelingAgent',
    'FigureFactoryAgent',
    'ScientificWriterAgent',
    'ReferenceAgent',
    'ReviewerAgent',
    'SubmissionAgent',
    'get_agent',
    'list_agents',
    'AGENT_REGISTRY'
]

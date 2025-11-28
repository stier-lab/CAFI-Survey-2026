"""
Project State Management
"""
from typing import Dict, List, Any, Optional
from datetime import datetime
from enum import Enum
import json
import os


class ProjectStage(Enum):
    """Project workflow stages."""
    INITIALIZATION = "initialization"
    FRAMING = "framing"
    LITERATURE = "literature"
    FRAMEWORK = "framework"
    DATA_QA = "data_qa"
    EDA = "eda"
    MODELING = "modeling"
    FIGURES = "figures"
    WRITING = "writing"
    REVIEW = "review"
    SUBMISSION = "submission"
    COMPLETED = "completed"


class TaskStatus(Enum):
    """Task completion status."""
    NOT_STARTED = "not_started"
    IN_PROGRESS = "in_progress"
    COMPLETED = "completed"
    BLOCKED = "blocked"


class ProjectState:
    """Manages the state of a research project."""

    def __init__(self, project_id: str, project_title: str):
        self.project_id = project_id
        self.project_title = project_title
        self.created_at = datetime.now().isoformat()
        self.updated_at = datetime.now().isoformat()
        self.current_stage = ProjectStage.INITIALIZATION
        self.project_brief: Optional[str] = None
        self.artifacts: Dict[str, Any] = {}
        self.tasks: List[Dict[str, Any]] = []
        self.agent_outputs: Dict[str, List[Dict[str, Any]]] = {}

    def set_project_brief(self, brief: str):
        """Set the project brief."""
        self.project_brief = brief
        self.updated_at = datetime.now().isoformat()

    def add_artifact(self, artifact_type: str, artifact_name: str,
                    content: Any, metadata: Optional[Dict[str, Any]] = None):
        """Add an artifact to the project."""
        if artifact_type not in self.artifacts:
            self.artifacts[artifact_type] = []

        self.artifacts[artifact_type].append({
            'name': artifact_name,
            'content': content,
            'metadata': metadata or {},
            'created_at': datetime.now().isoformat()
        })
        self.updated_at = datetime.now().isoformat()

    def get_artifacts(self, artifact_type: Optional[str] = None) -> List[Dict[str, Any]]:
        """Get artifacts, optionally filtered by type."""
        if artifact_type:
            return self.artifacts.get(artifact_type, [])
        return [item for items in self.artifacts.values() for item in items]

    def add_task(self, task_name: str, agent: str, description: str,
                dependencies: Optional[List[str]] = None):
        """Add a task to the project."""
        self.tasks.append({
            'id': f"task_{len(self.tasks) + 1}",
            'name': task_name,
            'agent': agent,
            'description': description,
            'status': TaskStatus.NOT_STARTED.value,
            'dependencies': dependencies or [],
            'created_at': datetime.now().isoformat(),
            'completed_at': None
        })
        self.updated_at = datetime.now().isoformat()

    def update_task_status(self, task_id: str, status: TaskStatus):
        """Update task status."""
        for task in self.tasks:
            if task['id'] == task_id:
                task['status'] = status.value
                if status == TaskStatus.COMPLETED:
                    task['completed_at'] = datetime.now().isoformat()
                self.updated_at = datetime.now().isoformat()
                break

    def add_agent_output(self, agent_name: str, output: Dict[str, Any]):
        """Store output from an agent execution."""
        if agent_name not in self.agent_outputs:
            self.agent_outputs[agent_name] = []

        self.agent_outputs[agent_name].append({
            'timestamp': datetime.now().isoformat(),
            'output': output
        })
        self.updated_at = datetime.now().isoformat()

    def get_agent_outputs(self, agent_name: str) -> List[Dict[str, Any]]:
        """Get all outputs from a specific agent."""
        return self.agent_outputs.get(agent_name, [])

    def set_stage(self, stage: ProjectStage):
        """Update project stage."""
        self.current_stage = stage
        self.updated_at = datetime.now().isoformat()

    def get_progress(self) -> Dict[str, Any]:
        """Calculate project progress."""
        total_tasks = len(self.tasks)
        if total_tasks == 0:
            return {'completion_percentage': 0, 'tasks_completed': 0, 'tasks_total': 0}

        completed_tasks = sum(1 for task in self.tasks
                            if task['status'] == TaskStatus.COMPLETED.value)

        return {
            'completion_percentage': (completed_tasks / total_tasks) * 100,
            'tasks_completed': completed_tasks,
            'tasks_total': total_tasks,
            'current_stage': self.current_stage.value
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize project state to dictionary."""
        return {
            'project_id': self.project_id,
            'project_title': self.project_title,
            'created_at': self.created_at,
            'updated_at': self.updated_at,
            'current_stage': self.current_stage.value,
            'project_brief': self.project_brief,
            'artifacts': self.artifacts,
            'tasks': self.tasks,
            'agent_outputs': self.agent_outputs
        }

    def save(self, directory: str):
        """Save project state to disk."""
        os.makedirs(directory, exist_ok=True)
        filepath = os.path.join(directory, f"{self.project_id}.json")

        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, filepath: str) -> 'ProjectState':
        """Load project state from disk."""
        with open(filepath, 'r') as f:
            data = json.load(f)

        project = cls(data['project_id'], data['project_title'])
        project.created_at = data['created_at']
        project.updated_at = data['updated_at']
        project.current_stage = ProjectStage(data['current_stage'])
        project.project_brief = data['project_brief']
        project.artifacts = data['artifacts']
        project.tasks = data['tasks']
        project.agent_outputs = data['agent_outputs']

        return project

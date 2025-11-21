"""
Orchestrator Agent - Coordinates all specialized agents
"""
from typing import Dict, List, Any
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.base_agent import BaseAgent
from core.project_state import ProjectState, ProjectStage, TaskStatus


class OrchestratorAgent(BaseAgent):
    """Orchestrates the entire research workflow."""

    def __init__(self):
        super().__init__(
            name="Orchestrator Agent",
            role="Coordinate all specialized agents to move a research project from idea to submission-ready manuscript",
            mission="Plan and oversee the sequence of steps across all agents, track artifacts, and ensure coherence between question, data, analyses, figures, and writing"
        )

    def get_system_prompt(self) -> str:
        return """You are the **Orchestrator Agent** for an academic research project.

Your job is to coordinate a suite of specialized agents to take a project from idea to submission-ready manuscript.

The available agents are:
1. Research PRD Agent - Converts research question into structured problem definition
2. Deep Literature Agent - Conducts literature exploration and synthesis
3. Conceptual Framework Agent - Creates conceptual diagrams and frameworks
4. Data QA Agent - Designs data quality assessment plan
5. EDA Agent - Explores dataset systematically
6. Modeling Agent - Designs statistical analysis plan
7. Figure Factory Agent - Designs all figures and visualizations
8. Scientific Writer Agent - Writes manuscript sections
9. Reference Agent - Manages citations and references
10. Reviewer Agent - Provides internal review feedback
11. Submission Agent - Prepares submission materials

Always:
- Start from the "Research Project Brief"
- Identify the current stage of work
- Decide which specialized agents should act next
- Write clear, concise task instructions for those agents including:
  * What they should read
  * What they should produce
  * How their output will be used in later stages

Return:
1. A short summary of the current project state
2. A list of 1-3 next actions with the agent name, input, and expected output
3. A project checklist with completed / in-progress / not-started tasks

Do not rewrite specialized outputs in detail (leave that to the relevant agent).
Do not invent data or results.
Always preserve traceability between decisions and original prompts."""

    def get_inputs(self) -> List[str]:
        return ["project_brief", "current_state"]

    def get_outputs(self) -> List[str]:
        return ["project_summary", "next_actions", "project_checklist"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Analyze project state and determine next steps.
        """
        try:
            is_valid, error_msg = self.validate_inputs(inputs)
            if not is_valid:
                raise ValueError(error_msg)

            project_brief = inputs.get('project_brief', '')
            current_state = inputs.get('current_state', {})

            # Analyze current stage and determine next steps
            next_actions = self._plan_next_actions(project_brief, current_state)

            # Generate project summary
            summary = self._generate_summary(current_state)

            # Create checklist
            checklist = self._generate_checklist(current_state)

            outputs = {
                'project_summary': summary,
                'next_actions': next_actions,
                'project_checklist': checklist
            }

            self.log_execution(inputs, outputs, 'success')
            return outputs

        except Exception as e:
            self.log_execution(inputs, {}, 'error', str(e))
            raise

    def _plan_next_actions(self, project_brief: str, current_state: Dict[str, Any]) -> List[Dict[str, str]]:
        """Determine the next 1-3 actions based on current state."""
        stage = current_state.get('current_stage', 'initialization')
        completed_agents = current_state.get('completed_agents', [])

        actions = []

        # Decision tree for workflow
        if stage == 'initialization' or 'research_prd' not in completed_agents:
            actions.append({
                'agent': 'Research PRD Agent',
                'task': 'Create Research Product Requirements Document',
                'inputs': ['project_brief'],
                'outputs': ['research_prd'],
                'rationale': 'Need to convert high-level question into precise, structured problem definition'
            })

        elif 'deep_literature' not in completed_agents:
            actions.append({
                'agent': 'Deep Literature Agent',
                'task': 'Conduct structured literature exploration',
                'inputs': ['project_brief', 'research_prd'],
                'outputs': ['literature_synthesis', 'gap_analysis', 'intro_outline'],
                'rationale': 'Need to anchor question in literature and identify gaps'
            })

        elif 'conceptual_framework' not in completed_agents:
            actions.append({
                'agent': 'Conceptual Framework Agent',
                'task': 'Create conceptual frameworks and diagrams',
                'inputs': ['research_prd', 'literature_synthesis'],
                'outputs': ['conceptual_frameworks', 'figure1_design'],
                'rationale': 'Need to create logical bridge between literature and hypotheses'
            })

        elif 'data_qa' not in completed_agents:
            actions.append({
                'agent': 'Data QA Agent',
                'task': 'Design data quality assessment plan',
                'inputs': ['research_prd', 'data_dictionary'],
                'outputs': ['qa_playbook', 'qa_checks', 'qa_report_template'],
                'rationale': 'Need to ensure data is structurally sound before analysis'
            })

        elif 'eda' not in completed_agents:
            actions.append({
                'agent': 'EDA Agent',
                'task': 'Systematically explore the dataset',
                'inputs': ['research_prd', 'cleaned_data', 'qa_report'],
                'outputs': ['eda_plan', 'eda_plots', 'eda_interpretations'],
                'rationale': 'Need to characterize distributions and patterns to inform modeling'
            })

        elif 'modeling' not in completed_agents:
            actions.append({
                'agent': 'Modeling Agent',
                'task': 'Design statistical analysis plan',
                'inputs': ['research_prd', 'eda_outputs'],
                'outputs': ['modeling_prd', 'model_specifications', 'code_skeletons'],
                'rationale': 'Need defensible modeling strategy aligned with hypotheses'
            })

        elif 'figures' not in completed_agents:
            actions.append({
                'agent': 'Figure Factory Agent',
                'task': 'Design all main and supplementary figures',
                'inputs': ['research_prd', 'modeling_plan', 'eda_insights'],
                'outputs': ['figure_plan', 'panel_designs', 'caption_drafts'],
                'rationale': 'Need coherent figure plan tied to hypotheses'
            })

        elif 'writing' not in completed_agents:
            # Can parallelize writing different sections
            actions.extend([
                {
                    'agent': 'Scientific Writer Agent',
                    'task': 'Write Introduction section',
                    'inputs': ['literature_synthesis', 'conceptual_framework'],
                    'outputs': ['draft_introduction'],
                    'rationale': 'Need to motivate the question and hypotheses'
                },
                {
                    'agent': 'Scientific Writer Agent',
                    'task': 'Write Methods section',
                    'inputs': ['modeling_plan', 'qa_playbook'],
                    'outputs': ['draft_methods'],
                    'rationale': 'Need to describe study design and analyses'
                }
            ])

        elif 'review' not in completed_agents:
            actions.append({
                'agent': 'Reviewer Agent',
                'task': 'Conduct internal review of manuscript',
                'inputs': ['full_manuscript', 'figures', 'modeling_plan'],
                'outputs': ['reviewer_reports', 'revision_suggestions'],
                'rationale': 'Need to identify weaknesses before external submission'
            })

        elif 'submission' not in completed_agents:
            actions.append({
                'agent': 'Submission Agent',
                'task': 'Prepare submission materials',
                'inputs': ['final_manuscript', 'figures', 'journal_guidelines'],
                'outputs': ['cover_letter', 'highlights', 'author_statements'],
                'rationale': 'Need to package manuscript for journal submission'
            })

        return actions

    def _generate_summary(self, current_state: Dict[str, Any]) -> str:
        """Generate a summary of current project state."""
        stage = current_state.get('current_stage', 'initialization')
        progress = current_state.get('progress', {})
        completed_agents = current_state.get('completed_agents', [])

        summary = f"**Current Stage:** {stage.replace('_', ' ').title()}\n\n"
        summary += f"**Progress:** {progress.get('completion_percentage', 0):.1f}% complete "
        summary += f"({progress.get('tasks_completed', 0)}/{progress.get('tasks_total', 0)} tasks)\n\n"
        summary += f"**Completed Agents:** {', '.join(completed_agents) if completed_agents else 'None yet'}\n"

        return summary

    def _generate_checklist(self, current_state: Dict[str, Any]) -> List[Dict[str, str]]:
        """Generate project checklist."""
        completed_agents = current_state.get('completed_agents', [])

        checklist_items = [
            {'task': 'Research PRD', 'status': 'completed' if 'research_prd' in completed_agents else 'not_started'},
            {'task': 'Literature Synthesis', 'status': 'completed' if 'deep_literature' in completed_agents else 'not_started'},
            {'task': 'Conceptual Framework', 'status': 'completed' if 'conceptual_framework' in completed_agents else 'not_started'},
            {'task': 'Data QA', 'status': 'completed' if 'data_qa' in completed_agents else 'not_started'},
            {'task': 'EDA', 'status': 'completed' if 'eda' in completed_agents else 'not_started'},
            {'task': 'Modeling', 'status': 'completed' if 'modeling' in completed_agents else 'not_started'},
            {'task': 'Figures', 'status': 'completed' if 'figures' in completed_agents else 'not_started'},
            {'task': 'Writing', 'status': 'completed' if 'writing' in completed_agents else 'not_started'},
            {'task': 'Review', 'status': 'completed' if 'review' in completed_agents else 'not_started'},
            {'task': 'Submission', 'status': 'completed' if 'submission' in completed_agents else 'not_started'},
        ]

        return checklist_items


if __name__ == "__main__":
    # Test the orchestrator
    orchestrator = OrchestratorAgent()
    print(f"Agent: {orchestrator.name}")
    print(f"Role: {orchestrator.role}")
    print(f"Inputs: {orchestrator.get_inputs()}")
    print(f"Outputs: {orchestrator.get_outputs()}")

"""
Research Workflow Dashboard - Streamlit Application
"""
import streamlit as st
import sys
import os
from datetime import datetime
from typing import Dict, Any, Optional

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.project_state import ProjectState, ProjectStage, TaskStatus
from agents.all_agents import get_agent, list_agents, AGENT_REGISTRY


# Page configuration
st.set_page_config(
    page_title="Research Workflow Dashboard",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)


def init_session_state():
    """Initialize session state variables."""
    if 'current_project' not in st.session_state:
        st.session_state.current_project = None
    if 'projects' not in st.session_state:
        st.session_state.projects = {}
    if 'selected_agent' not in st.session_state:
        st.session_state.selected_agent = None


def create_new_project():
    """Create a new research project."""
    st.header("Create New Research Project")

    with st.form("new_project_form"):
        project_title = st.text_input("Project Title", placeholder="e.g., Impact of Climate Change on Marine Biodiversity")
        project_id = st.text_input("Project ID", placeholder="e.g., marine-biodiv-2024")

        st.subheader("Project Brief")
        st.markdown("Fill out the research project brief that will guide all agents:")

        umbrella_question = st.text_area(
            "Umbrella Research Question",
            placeholder="What is the high-level research question?",
            height=100
        )

        focal_system = st.text_input("Focal System", placeholder="e.g., Coral reef ecosystems")
        key_taxa = st.text_input("Key Taxa/Entities", placeholder="e.g., Coral species, Fish communities")

        data_source = st.text_area(
            "Data Source(s)",
            placeholder="Describe your data sources",
            height=100
        )

        target_journals = st.text_input("Target Journals", placeholder="e.g., Nature Ecology & Evolution, Ecology Letters")

        submit_button = st.form_submit_button("Create Project")

        if submit_button:
            if not project_title or not project_id:
                st.error("Please provide both project title and ID")
                return

            # Create project brief
            project_brief = f"""# Research Project Brief

## 1. Project Title
{project_title}

## 2. Umbrella Question
{umbrella_question}

## 3. Nested Questions & Hypotheses
- Q1: [To be defined]
  - H1a: [To be defined]
  - H1b: [To be defined]
- Q2: [To be defined]

## 4. System & Concepts
- Focal system: {focal_system}
- Key taxa / entities: {key_taxa}
- Scales (space, time): [To be defined]
- Core mechanisms of interest: [To be defined]

## 5. Data Overview
- Data source(s): {data_source}
- Observational vs experimental: [To be defined]
- Response variables: [To be defined]
- Predictor variables: [To be defined]
- Random effects / grouping: [To be defined]

## 6. Target Outlet & Constraints
- Target journals: {target_journals}
- Word limits: [Check journal guidelines]
- Figure/table limits: [Check journal guidelines]

## 7. Deliverables
- Main manuscript
- Supplementary information
- Core figures (F1–F4 or more)
- Supplementary figures/tables
- Reproducible code repo
- Data & code availability statements

## 8. Known Risks & Confounders
[To be identified]
"""

            # Create new project
            project = ProjectState(project_id, project_title)
            project.set_project_brief(project_brief)

            # Add initial tasks
            project.add_task(
                "Create Research PRD",
                "research_prd",
                "Convert project brief into structured Research PRD",
                []
            )

            # Save project
            st.session_state.projects[project_id] = project
            st.session_state.current_project = project_id

            # Save to disk
            os.makedirs("data/projects", exist_ok=True)
            project.save("data/projects")

            st.success(f"Project '{project_title}' created successfully!")
            st.rerun()


def display_project_overview():
    """Display overview of current project."""
    if not st.session_state.current_project:
        st.info("No project selected. Create a new project or load an existing one.")
        return

    project = st.session_state.projects[st.session_state.current_project]

    st.header(f"📊 {project.project_title}")

    # Project metadata
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Project ID", project.project_id)
    with col2:
        st.metric("Current Stage", project.current_stage.value.replace('_', ' ').title())
    with col3:
        progress = project.get_progress()
        st.metric("Progress", f"{progress['completion_percentage']:.1f}%")

    # Progress bar
    st.progress(progress['completion_percentage'] / 100)
    st.caption(f"{progress['tasks_completed']} of {progress['tasks_total']} tasks completed")

    # Project brief
    with st.expander("View Project Brief", expanded=False):
        st.markdown(project.project_brief)

    # Task checklist
    st.subheader("Task Checklist")
    tasks_df_data = []
    for task in project.tasks:
        status_emoji = {
            'not_started': '⚪',
            'in_progress': '🟡',
            'completed': '✅',
            'blocked': '🔴'
        }
        tasks_df_data.append({
            'Status': status_emoji.get(task['status'], '⚪'),
            'Task': task['name'],
            'Agent': task['agent'],
            'Description': task['description']
        })

    if tasks_df_data:
        st.table(tasks_df_data)
    else:
        st.info("No tasks created yet. Use the Orchestrator to plan your workflow.")


def agent_interaction_panel():
    """Panel for interacting with agents."""
    if not st.session_state.current_project:
        st.warning("Please create or select a project first.")
        return

    st.header("🤖 Agent Interaction")

    project = st.session_state.projects[st.session_state.current_project]

    # Agent selector
    agents = list_agents()
    agent_names = [agent['name'] for agent in agents]
    agent_types = {agent['name']: agent['type'] for agent in agents}

    # Check if we should preselect an agent (from "Go to Agent" button)
    default_index = 0
    if 'preselect_agent' in st.session_state and st.session_state.preselect_agent:
        preselect_name = st.session_state.preselect_agent
        # Find matching agent name
        for i, name in enumerate(agent_names):
            if preselect_name in name or name in preselect_name:
                default_index = i
                break
        # Clear the preselect after using it
        st.session_state.preselect_agent = None

    selected_agent_name = st.selectbox(
        "Select Agent",
        options=agent_names,
        index=default_index,
        help="Choose which agent to work with"
    )

    agent_type = agent_types[selected_agent_name]
    agent = get_agent(agent_type)

    # Display agent info
    st.subheader(f"{agent.name}")
    st.markdown(f"**Role:** {agent.role}")
    st.markdown(f"**Mission:** {agent.mission}")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Required Inputs:**")
        for inp in agent.get_inputs():
            st.markdown(f"- {inp}")

    with col2:
        st.markdown("**Expected Outputs:**")
        for out in agent.get_outputs():
            st.markdown(f"- {out}")

    # System prompt
    with st.expander("View System Prompt", expanded=False):
        st.code(agent.get_system_prompt(), language="markdown")

    # Input collection
    st.subheader("Provide Inputs")
    inputs = {}

    for input_name in agent.get_inputs():
        # Check if input is available from previous agents
        available_artifacts = project.get_artifacts(input_name)

        if input_name == 'project_brief':
            inputs[input_name] = project.project_brief
            st.success(f"✓ {input_name} loaded from project")
        elif input_name == 'current_state':
            # Special handling for current_state - build from project
            current_state = {
                'current_stage': project.current_stage.value,
                'completed_agents': list(project.agent_outputs.keys()),
                'progress': project.get_progress()
            }
            inputs[input_name] = current_state
            st.success(f"✓ {input_name} built from project state")
            with st.expander("View current_state"):
                st.json(current_state)
        elif available_artifacts:
            st.success(f"✓ {input_name} available from previous work")
            inputs[input_name] = available_artifacts[-1]['content']
        else:
            user_input = st.text_area(
                f"{input_name}",
                placeholder=f"Provide {input_name} or it will be marked as [Not Provided]",
                height=150,
                key=f"input_{agent_type}_{input_name}"
            )
            inputs[input_name] = user_input if user_input else f"[{input_name} not provided]"

    # Execute agent
    if st.button(f"Execute {agent.name}", type="primary"):
        # Create containers for progress feedback
        progress_container = st.container()
        results_container = st.container()

        with progress_container:
            progress_bar = st.progress(0)
            status_text = st.empty()

            try:
                # Step 1: Initialize
                status_text.info(f"🔄 Initializing {agent.name}...")
                progress_bar.progress(10)

                # Step 2: Prepare inputs
                status_text.info(f"📋 Preparing inputs...")
                progress_bar.progress(20)
                context = {'project_brief': project.project_brief}

                # Step 3: Execute agent (this is the long operation)
                status_text.info(f"⚙️ Executing {agent.name}... (this may take 30-60 seconds)")
                progress_bar.progress(30)
                outputs = agent.execute(inputs, context)
                progress_bar.progress(70)

                # Step 4: Store outputs
                status_text.info(f"💾 Saving artifacts...")
                progress_bar.progress(80)
                for output_name, output_content in outputs.items():
                    project.add_artifact(
                        output_name,
                        f"{output_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
                        output_content,
                        {'agent': agent.name, 'agent_type': agent_type}
                    )

                # Store agent output
                project.add_agent_output(agent_type, outputs)

                # Step 5: Save project
                status_text.info(f"💾 Saving project state...")
                progress_bar.progress(90)
                project.save("data/projects")

                # Complete
                progress_bar.progress(100)
                status_text.success(f"✅ {agent.name} executed successfully!")

                # Display outputs in results container
                with results_container:
                    st.subheader("Outputs")
                    for output_name, output_content in outputs.items():
                        with st.expander(f"📄 {output_name}", expanded=True):
                            # Handle both string and dict/object outputs
                            if isinstance(output_content, str):
                                st.markdown(output_content)
                            elif isinstance(output_content, (dict, list)):
                                st.json(output_content)
                            else:
                                st.write(output_content)

            except ValueError as e:
                progress_bar.empty()
                status_text.empty()
                st.error("❌ Input Validation Error")
                st.markdown(f"""
                **Issue:** {str(e)}

                **How to fix:**
                - Ensure all required inputs are provided
                - Check that input format matches expectations
                """)

            except TimeoutError as e:
                progress_bar.empty()
                status_text.empty()
                st.error("⏱️ Operation Timed Out")
                st.markdown("""
                The agent took too long to respond. Please try again.
                """)

            except Exception as e:
                progress_bar.empty()
                status_text.empty()
                st.error(f"❌ Error executing agent")
                with st.expander("Technical Details"):
                    st.code(f"{type(e).__name__}: {str(e)}")
                    st.exception(e)


def artifacts_viewer():
    """View all artifacts generated in the project."""
    if not st.session_state.current_project:
        st.warning("Please create or select a project first.")
        return

    st.header("📁 Project Artifacts")

    project = st.session_state.projects[st.session_state.current_project]
    artifacts = project.get_artifacts()

    if not artifacts:
        st.info("No artifacts generated yet. Execute agents to create artifacts.")
        return

    # Group artifacts by type
    artifact_types = {}
    for artifact_list in project.artifacts.values():
        for artifact in artifact_list:
            artifact_type = artifact['name'].rsplit('_', 2)[0] if '_' in artifact['name'] else 'other'
            if artifact_type not in artifact_types:
                artifact_types[artifact_type] = []
            artifact_types[artifact_type].append(artifact)

    # Display artifacts
    for artifact_type, artifacts_of_type in artifact_types.items():
        st.subheader(artifact_type.replace('_', ' ').title())

        for artifact in artifacts_of_type:
            with st.expander(f"{artifact['name']} - {artifact['created_at'][:19]}"):
                # Handle both string and dict/object content
                content = artifact['content']
                if isinstance(content, str):
                    st.markdown(content)
                elif isinstance(content, (dict, list)):
                    st.json(content)
                else:
                    st.write(content)

                if artifact.get('metadata'):
                    st.caption(f"Generated by: {artifact['metadata'].get('agent', 'Unknown')}")

                # Download button - handle different content types
                download_content = content
                file_ext = "md"
                mime_type = "text/markdown"

                if isinstance(content, (dict, list)):
                    import json
                    download_content = json.dumps(content, indent=2)
                    file_ext = "json"
                    mime_type = "application/json"
                elif not isinstance(content, str):
                    download_content = str(content)

                st.download_button(
                    label="Download",
                    data=download_content,
                    file_name=f"{artifact['name']}.{file_ext}",
                    mime=mime_type,
                    key=f"download_{artifact['name']}_{artifact['created_at']}"
                )


def orchestrator_panel():
    """Special panel for orchestrator agent."""
    if not st.session_state.current_project:
        st.warning("Please create or select a project first.")
        return

    st.header("🎯 Orchestrator - Workflow Planning")

    project = st.session_state.projects[st.session_state.current_project]

    st.markdown("""
    The Orchestrator analyzes your project state and recommends the next steps.
    It coordinates all specialized agents to move your research from idea to submission.
    """)

    if st.button("Plan Next Steps", type="primary"):
        with st.spinner("Orchestrator analyzing project state..."):
            try:
                orchestrator = get_agent('orchestrator')

                # Prepare current state
                current_state = {
                    'current_stage': project.current_stage.value,
                    'completed_agents': list(project.agent_outputs.keys()),
                    'progress': project.get_progress()
                }

                inputs = {
                    'project_brief': project.project_brief,
                    'current_state': current_state
                }

                context = {'project_brief': project.project_brief}
                outputs = orchestrator.execute(inputs, context)

                # Save outputs to session state so they persist across reruns
                st.session_state.orchestrator_outputs = outputs

            except Exception as e:
                st.error(f"Error running orchestrator: {str(e)}")
                st.exception(e)

    # Display orchestrator outputs if available (outside the button block)
    if 'orchestrator_outputs' in st.session_state and st.session_state.orchestrator_outputs:
        outputs = st.session_state.orchestrator_outputs

        st.subheader("Project Summary")
        st.markdown(outputs['project_summary'])

        st.subheader("Recommended Next Actions")
        for i, action in enumerate(outputs['next_actions'], 1):
            with st.expander(f"Action {i}: {action['task']}", expanded=True):
                st.markdown(f"**Agent:** {action['agent']}")
                st.markdown(f"**Rationale:** {action['rationale']}")
                st.markdown(f"**Inputs Required:** {', '.join(action['inputs'])}")
                st.markdown(f"**Expected Outputs:** {', '.join(action['outputs'])}")

                if st.button(f"Go to {action['agent']}", key=f"goto_{i}"):
                    # Store the agent name to pre-select in agent interaction
                    st.session_state.preselect_agent = action['agent']
                    # Navigate to agent interaction page
                    st.session_state.page = "agents"
                    st.rerun()

        st.subheader("Project Checklist")
        checklist_df = []
        for item in outputs['project_checklist']:
            status_emoji = {
                'completed': '✅',
                'in_progress': '🟡',
                'not_started': '⚪'
            }
            checklist_df.append({
                'Status': status_emoji.get(item['status'], '⚪'),
                'Task': item['task']
            })
        st.table(checklist_df)


def sidebar_navigation():
    """Sidebar navigation."""
    with st.sidebar:
        st.title("🔬 Research Workflow")

        # Project selector
        st.subheader("Current Project")

        if st.session_state.projects:
            project_options = {
                project.project_title: project_id
                for project_id, project in st.session_state.projects.items()
            }

            if st.session_state.current_project:
                current_title = st.session_state.projects[st.session_state.current_project].project_title
                default_index = list(project_options.keys()).index(current_title)
            else:
                default_index = 0

            selected_project = st.selectbox(
                "Select Project",
                options=list(project_options.keys()),
                index=default_index
            )

            st.session_state.current_project = project_options[selected_project]

        if st.button("➕ New Project", use_container_width=True):
            st.session_state.page = "new_project"
            st.rerun()

        st.divider()

        # Navigation
        st.subheader("Navigation")

        pages = {
            "Overview": "overview",
            "Orchestrator": "orchestrator",
            "Agent Interaction": "agents",
            "Artifacts": "artifacts",
        }

        for page_name, page_id in pages.items():
            if st.button(page_name, use_container_width=True):
                st.session_state.page = page_id
                st.rerun()

        st.divider()

        # Quick agent access
        if st.session_state.current_project:
            st.subheader("Quick Agent Access")
            agents = list_agents()
            for agent in agents[:5]:  # Show first 5 agents
                if st.button(f"🤖 {agent['name']}", use_container_width=True, key=f"quick_{agent['type']}"):
                    st.session_state.page = "agents"
                    st.session_state.preselect_agent = agent['name']
                    st.rerun()


def main():
    """Main application."""
    init_session_state()

    # Initialize page state
    if 'page' not in st.session_state:
        st.session_state.page = "overview"

    # Load existing projects
    if os.path.exists("data/projects"):
        for filename in os.listdir("data/projects"):
            if filename.endswith(".json"):
                filepath = os.path.join("data/projects", filename)
                try:
                    project = ProjectState.load(filepath)
                    if project.project_id not in st.session_state.projects:
                        st.session_state.projects[project.project_id] = project
                except Exception as e:
                    st.error(f"Error loading project {filename}: {e}")

    # Sidebar
    sidebar_navigation()

    # Main content
    page = st.session_state.page

    if page == "new_project":
        create_new_project()
    elif page == "overview":
        display_project_overview()
    elif page == "orchestrator":
        orchestrator_panel()
    elif page == "agents":
        agent_interaction_panel()
    elif page == "artifacts":
        artifacts_viewer()
    else:
        st.error("Unknown page")


if __name__ == "__main__":
    main()

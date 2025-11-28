"""
Research Workflow Dashboard - Redesigned with UX Best Practices
Following Frontend Master Prompt v2.0 principles
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


# ============================================================================
# DESIGN SYSTEM CONFIGURATION
# ============================================================================

# Page configuration with proper meta
st.set_page_config(
    page_title="Research Workflow System",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': 'https://github.com/research-workflow/docs',
        'Report a bug': 'https://github.com/research-workflow/issues',
        'About': 'AI-powered research workflow management system'
    }
)

# Design tokens and CSS variables
st.markdown("""
<style>
    /* Design System Tokens */
    :root {
        /* Color Palette */
        --color-primary: #2563eb;
        --color-primary-hover: #1d4ed8;
        --color-secondary: #64748b;
        --color-success: #10b981;
        --color-warning: #f59e0b;
        --color-error: #ef4444;
        --color-info: #3b82f6;

        /* Neutral Scale */
        --color-background: #ffffff;
        --color-surface: #f8fafc;
        --color-surface-elevated: #ffffff;
        --color-border: #e2e8f0;
        --color-border-strong: #cbd5e1;

        /* Text Colors */
        --color-text-primary: #0f172a;
        --color-text-secondary: #475569;
        --color-text-tertiary: #64748b;
        --color-text-inverse: #ffffff;

        /* Spacing Scale (4px base) */
        --space-1: 0.25rem;
        --space-2: 0.5rem;
        --space-3: 0.75rem;
        --space-4: 1rem;
        --space-6: 1.5rem;
        --space-8: 2rem;
        --space-12: 3rem;
        --space-16: 4rem;

        /* Border Radius */
        --radius-sm: 4px;
        --radius-md: 6px;
        --radius-lg: 8px;
        --radius-xl: 12px;

        /* Shadows */
        --shadow-sm: 0 1px 2px rgba(0,0,0,0.05);
        --shadow-md: 0 4px 6px rgba(0,0,0,0.07), 0 2px 4px rgba(0,0,0,0.06);
        --shadow-lg: 0 10px 15px rgba(0,0,0,0.1), 0 4px 6px rgba(0,0,0,0.05);

        /* Typography */
        --font-sans: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        --font-mono: ui-monospace, "Cascadia Code", "Source Code Pro", Menlo, monospace;

        /* Transitions */
        --transition-fast: 150ms ease-out;
        --transition-base: 200ms ease-out;
        --transition-slow: 300ms ease-out;
    }

    /* Dark Mode Support */
    @media (prefers-color-scheme: dark) {
        :root {
            --color-background: #0f172a;
            --color-surface: #1e293b;
            --color-surface-elevated: #334155;
            --color-border: #334155;
            --color-border-strong: #475569;
            --color-text-primary: #f1f5f9;
            --color-text-secondary: #cbd5e1;
            --color-text-tertiary: #94a3b8;
        }
    }

    /* Motion Accessibility */
    @media (prefers-reduced-motion: reduce) {
        * {
            animation-duration: 0.01ms !important;
            animation-iteration-count: 1 !important;
            transition-duration: 0.01ms !important;
        }
    }

    /* Global Styles */
    .main {
        background-color: var(--color-background);
        color: var(--color-text-primary);
        font-family: var(--font-sans);
    }

    /* Hide Streamlit Branding */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}

    /* Typography Improvements */
    h1 {
        font-size: 2.25rem;
        font-weight: 700;
        line-height: 1.25;
        color: var(--color-text-primary);
        margin-bottom: var(--space-6);
    }

    h2 {
        font-size: 1.875rem;
        font-weight: 600;
        line-height: 1.25;
        color: var(--color-text-primary);
        margin-bottom: var(--space-4);
    }

    h3 {
        font-size: 1.5rem;
        font-weight: 600;
        line-height: 1.375;
        color: var(--color-text-primary);
        margin-bottom: var(--space-3);
    }

    p {
        line-height: 1.5;
        color: var(--color-text-secondary);
        margin-bottom: var(--space-4);
    }

    /* Improved Sidebar */
    [data-testid="stSidebar"] {
        background-color: var(--color-surface);
        border-right: 1px solid var(--color-border);
        padding: var(--space-6);
    }

    [data-testid="stSidebar"] h1 {
        font-size: 1.25rem;
        font-weight: 700;
        margin-bottom: var(--space-6);
        color: var(--color-text-primary);
    }

    /* Card Component */
    .card {
        background: var(--color-surface-elevated);
        border: 1px solid var(--color-border);
        border-radius: var(--radius-lg);
        padding: var(--space-6);
        margin-bottom: var(--space-4);
        box-shadow: var(--shadow-sm);
        transition: all var(--transition-base);
    }

    .card:hover {
        box-shadow: var(--shadow-md);
        transform: translateY(-2px);
    }

    /* Button Improvements */
    .stButton > button {
        background-color: var(--color-primary);
        color: var(--color-text-inverse);
        border: none;
        border-radius: var(--radius-md);
        padding: var(--space-2) var(--space-4);
        font-weight: 500;
        font-size: 0.875rem;
        transition: all var(--transition-fast);
        cursor: pointer;
        box-shadow: var(--shadow-sm);
    }

    .stButton > button:hover {
        background-color: var(--color-primary-hover);
        box-shadow: var(--shadow-md);
        transform: translateY(-1px);
    }

    .stButton > button:active {
        transform: scale(0.98);
    }

    .stButton > button:focus-visible {
        outline: 2px solid var(--color-primary);
        outline-offset: 2px;
    }

    /* Secondary Button */
    .stButton.secondary > button {
        background-color: var(--color-surface);
        color: var(--color-text-primary);
        border: 1px solid var(--color-border);
    }

    .stButton.secondary > button:hover {
        background-color: var(--color-surface-elevated);
        border-color: var(--color-border-strong);
    }

    /* Input Fields */
    .stTextInput > div > div > input {
        border: 1px solid var(--color-border);
        border-radius: var(--radius-md);
        padding: var(--space-2) var(--space-3);
        font-size: 0.875rem;
        transition: all var(--transition-fast);
        background-color: var(--color-surface-elevated);
        color: var(--color-text-primary);
    }

    .stTextInput > div > div > input:focus {
        border-color: var(--color-primary);
        box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.1);
        outline: none;
    }

    /* Text Area */
    .stTextArea > div > div > textarea {
        border: 1px solid var(--color-border);
        border-radius: var(--radius-md);
        padding: var(--space-3);
        font-size: 0.875rem;
        transition: all var(--transition-fast);
        background-color: var(--color-surface-elevated);
        color: var(--color-text-primary);
        line-height: 1.5;
    }

    .stTextArea > div > div > textarea:focus {
        border-color: var(--color-primary);
        box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.1);
        outline: none;
    }

    /* Select Box */
    .stSelectbox > div > div {
        border: 1px solid var(--color-border);
        border-radius: var(--radius-md);
        background-color: var(--color-surface-elevated);
    }

    /* Metric Cards */
    [data-testid="stMetricValue"] {
        font-size: 2rem;
        font-weight: 700;
        color: var(--color-text-primary);
    }

    [data-testid="stMetricLabel"] {
        font-size: 0.875rem;
        font-weight: 500;
        color: var(--color-text-secondary);
        text-transform: uppercase;
        letter-spacing: 0.025em;
    }

    /* Progress Bar */
    .stProgress > div > div {
        background-color: var(--color-border);
        border-radius: var(--radius-lg);
        height: 8px;
    }

    .stProgress > div > div > div {
        background-color: var(--color-primary);
        border-radius: var(--radius-lg);
        transition: width var(--transition-slow);
    }

    /* Status Badges */
    .status-badge {
        display: inline-flex;
        align-items: center;
        padding: var(--space-1) var(--space-3);
        border-radius: var(--radius-md);
        font-size: 0.75rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.025em;
    }

    .status-success {
        background-color: rgba(16, 185, 129, 0.1);
        color: var(--color-success);
    }

    .status-warning {
        background-color: rgba(245, 158, 11, 0.1);
        color: var(--color-warning);
    }

    .status-error {
        background-color: rgba(239, 68, 68, 0.1);
        color: var(--color-error);
    }

    .status-info {
        background-color: rgba(59, 130, 246, 0.1);
        color: var(--color-info);
    }

    /* Expander */
    .streamlit-expanderHeader {
        background-color: var(--color-surface);
        border: 1px solid var(--color-border);
        border-radius: var(--radius-md);
        padding: var(--space-3);
        font-weight: 500;
        transition: all var(--transition-fast);
    }

    .streamlit-expanderHeader:hover {
        background-color: var(--color-surface-elevated);
        border-color: var(--color-border-strong);
    }

    /* Toast/Alert */
    .alert {
        padding: var(--space-4);
        border-radius: var(--radius-md);
        margin-bottom: var(--space-4);
        border-left: 4px solid;
    }

    .alert-info {
        background-color: rgba(59, 130, 246, 0.1);
        border-color: var(--color-info);
        color: var(--color-text-primary);
    }

    .alert-success {
        background-color: rgba(16, 185, 129, 0.1);
        border-color: var(--color-success);
        color: var(--color-text-primary);
    }

    .alert-warning {
        background-color: rgba(245, 158, 11, 0.1);
        border-color: var(--color-warning);
        color: var(--color-text-primary);
    }

    .alert-error {
        background-color: rgba(239, 68, 68, 0.1);
        border-color: var(--color-error);
        color: var(--color-text-primary);
    }

    /* Loading State */
    .loading-skeleton {
        background: linear-gradient(
            90deg,
            var(--color-surface) 0%,
            var(--color-surface-elevated) 50%,
            var(--color-surface) 100%
        );
        background-size: 200% 100%;
        animation: loading 1.5s ease-in-out infinite;
        border-radius: var(--radius-md);
        height: 20px;
        margin-bottom: var(--space-2);
    }

    @keyframes loading {
        0% { background-position: 200% 0; }
        100% { background-position: -200% 0; }
    }

    /* Empty State */
    .empty-state {
        text-align: center;
        padding: var(--space-12);
        color: var(--color-text-secondary);
    }

    .empty-state-icon {
        font-size: 3rem;
        margin-bottom: var(--space-4);
        opacity: 0.5;
    }

    .empty-state-title {
        font-size: 1.25rem;
        font-weight: 600;
        color: var(--color-text-primary);
        margin-bottom: var(--space-2);
    }

    .empty-state-description {
        font-size: 0.875rem;
        color: var(--color-text-secondary);
        margin-bottom: var(--space-6);
    }

    /* Table Improvements */
    table {
        width: 100%;
        border-collapse: separate;
        border-spacing: 0;
    }

    thead th {
        background-color: var(--color-surface);
        border-bottom: 2px solid var(--color-border-strong);
        padding: var(--space-3);
        text-align: left;
        font-weight: 600;
        font-size: 0.875rem;
        color: var(--color-text-primary);
        text-transform: uppercase;
        letter-spacing: 0.025em;
    }

    tbody td {
        padding: var(--space-3);
        border-bottom: 1px solid var(--color-border);
        color: var(--color-text-secondary);
    }

    tbody tr:hover {
        background-color: var(--color-surface);
        transition: background-color var(--transition-fast);
    }

    /* Accessibility: Focus Indicators */
    *:focus-visible {
        outline: 2px solid var(--color-primary);
        outline-offset: 2px;
        border-radius: var(--radius-sm);
    }

    /* Skip Link for Accessibility */
    .skip-link {
        position: absolute;
        top: -40px;
        left: 0;
        background: var(--color-primary);
        color: var(--color-text-inverse);
        padding: var(--space-2) var(--space-4);
        text-decoration: none;
        border-radius: var(--radius-md);
        z-index: 1000;
    }

    .skip-link:focus {
        top: var(--space-2);
    }
</style>
""", unsafe_allow_html=True)


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def init_session_state():
    """Initialize session state with proper defaults."""
    if 'current_project' not in st.session_state:
        st.session_state.current_project = None
    if 'projects' not in st.session_state:
        st.session_state.projects = {}
    if 'selected_agent' not in st.session_state:
        st.session_state.selected_agent = None
    if 'show_success' not in st.session_state:
        st.session_state.show_success = False
    if 'success_message' not in st.session_state:
        st.session_state.success_message = ''


def show_alert(message: str, alert_type: str = "info"):
    """Display an accessible alert message."""
    alert_classes = {
        "info": "alert alert-info",
        "success": "alert alert-success",
        "warning": "alert alert-warning",
        "error": "alert alert-error"
    }

    icons = {
        "info": "ℹ️",
        "success": "✓",
        "warning": "⚠️",
        "error": "✕"
    }

    st.markdown(f"""
        <div class="{alert_classes.get(alert_type, 'alert alert-info')}" role="alert">
            <strong>{icons.get(alert_type, "ℹ️")}</strong> {message}
        </div>
    """, unsafe_allow_html=True)


def show_loading_skeleton(count: int = 3):
    """Display loading skeleton for better perceived performance."""
    for _ in range(count):
        st.markdown('<div class="loading-skeleton"></div>', unsafe_allow_html=True)


def show_empty_state(icon: str, title: str, description: str, action_label: str = None, action_key: str = None):
    """Display an engaging empty state."""
    st.markdown(f"""
        <div class="empty-state">
            <div class="empty-state-icon">{icon}</div>
            <h3 class="empty-state-title">{title}</h3>
            <p class="empty-state-description">{description}</p>
        </div>
    """, unsafe_allow_html=True)

    if action_label and action_key:
        col1, col2, col3 = st.columns([1, 1, 1])
        with col2:
            return st.button(action_label, key=action_key, use_container_width=True)
    return False


def get_status_badge(status: str) -> str:
    """Return HTML for status badge."""
    status_map = {
        'completed': ('✓', 'status-success'),
        'in_progress': ('⟳', 'status-warning'),
        'not_started': ('○', 'status-info'),
        'blocked': ('✕', 'status-error')
    }

    icon, badge_class = status_map.get(status, ('○', 'status-info'))
    status_text = status.replace('_', ' ').title()

    return f'<span class="status-badge {badge_class}">{icon} {status_text}</span>'


# ============================================================================
# PAGE COMPONENTS
# ============================================================================

def create_new_project_page():
    """Improved new project creation form with better UX."""
    st.title("Create New Research Project")

    st.markdown("""
        <p style="font-size: 1.125rem; color: var(--color-text-secondary); margin-bottom: var(--space-8);">
            Let's set up your research project. We'll guide you through defining your research question,
            hypotheses, and data sources.
        </p>
    """, unsafe_allow_html=True)

    with st.form("new_project_form", clear_on_submit=False):
        # Step 1: Basic Information
        st.subheader("📋 Basic Information")

        col1, col2 = st.columns(2)
        with col1:
            project_title = st.text_input(
                "Project Title",
                placeholder="e.g., Impact of Climate Change on Marine Biodiversity",
                help="Give your project a descriptive title",
                max_chars=100
            )

        with col2:
            project_id = st.text_input(
                "Project ID",
                placeholder="e.g., marine-biodiv-2024",
                help="A unique identifier (lowercase, hyphens allowed)",
                max_chars=50
            )

        st.divider()

        # Step 2: Research Question
        st.subheader("🔬 Research Question")

        umbrella_question = st.text_area(
            "Umbrella Research Question",
            placeholder="What is the high-level research question you're investigating?",
            help="Your main research question that encompasses the entire study",
            height=100
        )

        st.divider()

        # Step 3: System & Data
        st.subheader("🌍 System & Data")

        col1, col2 = st.columns(2)
        with col1:
            focal_system = st.text_input(
                "Focal System",
                placeholder="e.g., Coral reef ecosystems",
                help="What system are you studying?"
            )

        with col2:
            key_taxa = st.text_input(
                "Key Taxa/Entities",
                placeholder="e.g., Coral species, Fish communities",
                help="Main organisms or entities in your study"
            )

        data_source = st.text_area(
            "Data Source(s)",
            placeholder="Describe your data sources, collection methods, and timeframe",
            help="Be as specific as possible about your data",
            height=120
        )

        st.divider()

        # Step 4: Publication Goals
        st.subheader("📝 Publication Goals")

        target_journals = st.text_input(
            "Target Journals",
            placeholder="e.g., Nature Ecology & Evolution, Ecology Letters",
            help="Where do you plan to submit this research?"
        )

        # Form submission
        st.divider()

        col1, col2, col3 = st.columns([1, 1, 1])
        with col2:
            submit_button = st.form_submit_button(
                "Create Project",
                use_container_width=True,
                type="primary"
            )

        if submit_button:
            # Validation
            if not project_title or not project_id:
                show_alert("Please provide both project title and ID", "error")
                return

            if not umbrella_question:
                show_alert("Please provide an umbrella research question", "warning")
                return

            # Create project brief
            project_brief = f"""# Research Project Brief

## 1. Project Title
{project_title}

## 2. Umbrella Question
{umbrella_question}

## 3. Nested Questions & Hypotheses
- Q1: [To be defined with Research PRD Agent]
  - H1a: [To be defined]
  - H1b: [To be defined]
- Q2: [To be defined]

## 4. System & Concepts
- Focal system: {focal_system or '[To be defined]'}
- Key taxa / entities: {key_taxa or '[To be defined]'}
- Scales (space, time): [To be defined]
- Core mechanisms of interest: [To be defined]

## 5. Data Overview
- Data source(s): {data_source or '[To be defined]'}
- Observational vs experimental: [To be defined]
- Response variables: [To be defined]
- Predictor variables: [To be defined]
- Random effects / grouping: [To be defined]

## 6. Target Outlet & Constraints
- Target journals: {target_journals or '[To be defined]'}
- Word limits: [Check journal guidelines]
- Figure/table limits: [Check journal guidelines]

## 7. Deliverables
- Main manuscript
- Supplementary information
- Core figures (F1–F4 or more)
- Reproducible code repo
- Data & code availability statements

## 8. Known Risks & Confounders
[To be identified with Research PRD Agent]
"""

            # Create project
            project = ProjectState(project_id, project_title)
            project.set_project_brief(project_brief)

            # Add initial task
            project.add_task(
                "Create Research PRD",
                "research_prd",
                "Convert project brief into structured Research PRD",
                []
            )

            # Save
            st.session_state.projects[project_id] = project
            st.session_state.current_project = project_id

            os.makedirs("data/projects", exist_ok=True)
            project.save("data/projects")

            # Success feedback
            st.session_state.show_success = True
            st.session_state.success_message = f"Project '{project_title}' created successfully!"

            st.rerun()


def display_project_overview():
    """Enhanced project overview with better visual hierarchy."""
    if not st.session_state.current_project:
        if show_empty_state(
            icon="📂",
            title="No Project Selected",
            description="Create a new research project or load an existing one to get started with your research workflow.",
            action_label="Create New Project",
            action_key="empty_create_project"
        ):
            st.session_state.page = "new_project"
            st.rerun()
        return

    project = st.session_state.projects[st.session_state.current_project]

    # Header with project title
    st.title(f"📊 {project.project_title}")

    # Success message if just created
    if st.session_state.show_success:
        show_alert(st.session_state.success_message, "success")
        st.session_state.show_success = False

    # Metrics row
    progress = project.get_progress()
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric(
            label="Project ID",
            value=project.project_id
        )

    with col2:
        st.metric(
            label="Current Stage",
            value=project.current_stage.value.replace('_', ' ').title()
        )

    with col3:
        st.metric(
            label="Progress",
            value=f"{progress['completion_percentage']:.0f}%"
        )

    with col4:
        st.metric(
            label="Tasks Completed",
            value=f"{progress['tasks_completed']}/{progress['tasks_total']}"
        )

    # Progress bar
    st.progress(progress['completion_percentage'] / 100)

    st.divider()

    # Two column layout
    col1, col2 = st.columns([2, 1])

    with col1:
        # Project brief
        st.subheader("Research Brief")
        with st.expander("View Full Project Brief", expanded=False):
            st.markdown(project.project_brief)

    with col2:
        # Task checklist
        st.subheader("Task Checklist")

        if project.tasks:
            for task in project.tasks:
                status_html = get_status_badge(task['status'])
                st.markdown(f"""
                    <div style="padding: var(--space-3); margin-bottom: var(--space-2); background: var(--color-surface); border-radius: var(--radius-md); border-left: 3px solid var(--color-primary);">
                        <div style="display: flex; justify-content: space-between; align-items: center;">
                            <strong>{task['name']}</strong>
                            {status_html}
                        </div>
                        <div style="font-size: 0.875rem; color: var(--color-text-secondary); margin-top: var(--space-1);">
                            {task['agent']}
                        </div>
                    </div>
                """, unsafe_allow_html=True)
        else:
            show_empty_state(
                icon="📝",
                title="No Tasks Yet",
                description="Use the Orchestrator to plan your workflow and generate tasks.",
                action_label=None,
                action_key=None
            )


def sidebar_navigation():
    """Enhanced sidebar with better organization and accessibility."""
    with st.sidebar:
        # Skip link for accessibility
        st.markdown('<a href="#main-content" class="skip-link">Skip to main content</a>', unsafe_allow_html=True)

        st.markdown("# 🔬 Research Workflow")

        # Project selector section
        st.markdown("### Current Project")

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
                index=default_index,
                label_visibility="collapsed"
            )

            st.session_state.current_project = project_options[selected_project]

        if st.button("➕ New Project", use_container_width=True, type="primary"):
            st.session_state.page = "new_project"
            st.rerun()

        st.divider()

        # Navigation section
        st.markdown("### Navigation")

        pages = {
            "Overview": ("overview", "📊"),
            "Orchestrator": ("orchestrator", "🎯"),
            "Agent Interaction": ("agents", "🤖"),
            "Artifacts": ("artifacts", "📁"),
        }

        for page_name, (page_id, icon) in pages.items():
            if st.button(f"{icon} {page_name}", use_container_width=True, key=f"nav_{page_id}"):
                st.session_state.page = page_id
                st.rerun()

        st.divider()

        # Quick agent access
        if st.session_state.current_project:
            st.markdown("### Quick Agent Access")
            agents = list_agents()
            for agent in agents[:5]:
                if st.button(
                    f"🤖 {agent['name']}",
                    use_container_width=True,
                    key=f"quick_{agent['type']}",
                    help=agent['role']
                ):
                    st.session_state.page = "agents"
                    st.session_state.selected_agent = agent['type']
                    st.rerun()


def main():
    """Main application with improved UX flow."""
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
                    show_alert(f"Error loading project {filename}: {e}", "error")

    # Sidebar
    sidebar_navigation()

    # Main content with proper ARIA landmark
    st.markdown('<div id="main-content" role="main">', unsafe_allow_html=True)

    # Route to appropriate page
    page = st.session_state.page

    if page == "new_project":
        create_new_project_page()
    elif page == "overview":
        display_project_overview()
    elif page == "orchestrator":
        st.title("🎯 Orchestrator - Workflow Planning")
        st.info("Orchestrator page coming soon with enhanced UX")
    elif page == "agents":
        st.title("🤖 Agent Interaction")
        st.info("Agent interaction page coming soon with enhanced UX")
    elif page == "artifacts":
        st.title("📁 Project Artifacts")
        st.info("Artifacts page coming soon with enhanced UX")
    else:
        show_alert("Unknown page", "error")

    st.markdown('</div>', unsafe_allow_html=True)


if __name__ == "__main__":
    main()

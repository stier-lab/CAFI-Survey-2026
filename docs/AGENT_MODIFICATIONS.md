# Agent System for CAFI Survey Analysis

This document describes the research workflow agent system and recommended modifications for the CAFI Survey Analysis project.

## System Overview

The agent system provides AI-assisted research workflow automation with 12 specialized agents covering the complete research lifecycle from problem definition to manuscript submission.

## Current Status: READY TO USE

All core modules have been copied from the research-agent project:
- `core/base_agent.py` - Abstract base class for agents
- `core/project_state.py` - Project state management
- `core/llm_interface.py` - LLM provider abstraction
- `config.py` - Configuration settings
- `dashboard/` - Streamlit web interface
- `run.py` - Main entry point

---

## Quick Start

### 1. Install Dependencies

```bash
cd /Users/adrianstiermbp2023/CAFI-Survey-2026
pip install -r requirements.txt
```

### 2. Configure API Key

Create a `.env` file or set environment variable:

```bash
export ANTHROPIC_API_KEY="your-api-key-here"
```

### 3. Run the Dashboard

```bash
python run.py
# Or directly:
streamlit run dashboard/app.py
```

### 4. Run Agents Programmatically

```python
from agents import get_agent
from core.project_state import ProjectState

# Initialize project
state = ProjectState(project_id="cafi-survey-2026", project_title="CAFI Survey Analysis")

# Run Research PRD Agent
prd_agent = get_agent('research_prd')
outputs = prd_agent.execute(
    inputs={
        'project_brief': open('docs/PRD.md').read(),
        'data_dictionary': 'See data/README files'
    },
    context={'use_llm': True}
)

print(outputs['research_prd'])
```

---

## Available Agents

| Agent | Purpose | Inputs | Outputs |
|-------|---------|--------|---------|
| **Orchestrator** | Workflow coordination | project_brief, current_state | next_actions, checklist |
| **Research PRD** | Problem definition | project_brief, data_dictionary | research_prd, hypotheses |
| **Literature** | Literature synthesis | research_prd, seed_papers | synthesis, gap_analysis |
| **Framework** | Conceptual diagrams | research_prd, literature | frameworks, figure_design |
| **Data QA** | Quality assessment | research_prd, data_dictionary | qa_playbook, checks |
| **EDA** | Exploratory analysis | research_prd, cleaned_data | eda_plan, plots |
| **Modeling** | Statistical planning | research_prd, eda_outputs | model_specs, code |
| **Figures** | Visualization design | research_prd, modeling_plan | figure_inventory, captions |
| **Writer** | Manuscript sections | all_outputs | introduction, methods, results |
| **References** | Citation management | manuscript, refs | formatted_refs |
| **Reviewer** | Internal review | manuscript, figures | feedback, revisions |
| **Submission** | Submission prep | final_manuscript | cover_letter, highlights |

---

## Recommended CAFI-Specific Modifications

### 1. Add Domain Context to Agent Prompts

Create a shared context string for all agents:

```python
CAFI_CONTEXT = """
## CAFI Project Context

**CAFI** = Crabs, Alpheid shrimp, Fish, and snails (Invertebrates)

**Study System:**
- Location: Mo'orea, French Polynesia
- Sites: HAU, MAT, MRB
- Host: Pocillopora spp.
- Sample: 114 corals, 2,847 CAFI, 10 OTUs

**Key Variables:**
- Response: abundance, richness, diversity
- Predictors: volume, depth, branch_width, morphotype
- Random effects: site

**Important Notes:**
- OTUs are morphological (no genetic confirmation)
- Branch width (tight/wide) is the functional trait
- Physiology data is subset only
"""
```

Then modify each agent's `get_system_prompt()` to include this context.

### 2. Create R Pipeline Bridge

Add a bridge module to read R outputs:

```python
# agents/r_bridge.py
import pandas as pd
from pathlib import Path

class RBridge:
    def __init__(self, project_root):
        self.root = Path(project_root)
        self.output = self.root / "output"

    def get_data_dictionary(self):
        path = self.output / "tables" / "data_dictionary.csv"
        if path.exists():
            return pd.read_csv(path).to_markdown()
        return "Run 01_load_clean_data.R first"

    def get_model_results(self):
        tables = self.output / "tables"
        return {
            'coefficients': list(tables.glob("*coefficients*.csv")),
            'importance': list(tables.glob("*importance*.csv"))
        }
```

### 3. Update DataQAAgent for CAFI Data

Modify the system prompt to include CAFI-specific checks:

```python
def get_system_prompt(self):
    return f"""{CAFI_CONTEXT}

You are the Data QA Agent for CAFI Survey Analysis.

Include checks for:
1. CAFI-coral merge (many-to-one on coral_id)
2. Variable ranges: volume 0-10000 cm³, depth 0-20m
3. Categorical levels: site (HAU/MAT/MRB), morphotype
4. Ecological sanity: positive size-abundance correlation
5. Missing data: physiology subset only
"""
```

### 4. Update ModelingAgent for Survey Statistics

```python
def get_system_prompt(self):
    return f"""{CAFI_CONTEXT}

You are the Modeling Agent for CAFI Survey Analysis.

Design models including:
1. GLMMs: abundance ~ volume + (1|site) [negative binomial]
2. PERMANOVA: community ~ site + morphotype
3. Spatial: Moran's I for autocorrelation
4. ML: Random Forest for variable importance

Address spatial autocorrelation in p-values.
Report effect sizes and 95% CIs.
"""
```

### 5. Update FigureFactoryAgent for Publication Standards

```python
def get_system_prompt(self):
    return f"""{CAFI_CONTEXT}

You are the Figure Factory Agent.

Standards:
- Resolution: 300 dpi, white background
- Colors: viridis palette
- Fonts: sans-serif, 10-12pt

Required figures:
1. Study overview (map, distributions)
2. Coral-CAFI relationships (scatterplots, boxplots)
3. Community composition (NMDS, accumulation)
4. Spatial patterns (Moran's I, LISA)

Return ggplot2 code outlines.
"""
```

---

## Integration with R Pipeline

The agents can be used alongside the R analysis scripts:

### Workflow
1. Run R scripts (00-17) to generate data outputs
2. Use agents to interpret and synthesize results
3. Agent outputs inform manuscript writing

### Data Flow
```
R Scripts → output/tables/*.csv → R Bridge → Agents
                                           ↓
                              → Manuscript sections
                              → Figure captions
                              → Statistical interpretations
```

### Example: Generate Methods Section

```python
from agents import get_agent
from agents.r_bridge import RBridge

bridge = RBridge('/path/to/CAFI-Survey-2026')
writer = get_agent('writer')

# Get model specifications from R outputs
model_results = bridge.get_model_results()

# Generate methods section
outputs = writer.execute(
    inputs={
        'modeling_plan': str(model_results),
        'data_dictionary': bridge.get_data_dictionary()
    },
    context={'requested_section': 'methods'}
)

print(outputs['draft_methods'])
```

---

## Project State Management

The system tracks project progress:

```python
from core.project_state import ProjectState, ProjectStage

# Initialize
state = ProjectState(
    project_id="cafi-survey-2026",
    project_title="CAFI Survey Analysis"
)

# Track artifacts
state.add_artifact('research_prd', 'docs/PRD.md')
state.add_artifact('data_dictionary', 'output/tables/data_dictionary.csv')

# Update progress
state.current_stage = ProjectStage.MODELING
state.add_agent_output('modeling', model_outputs)

# Save/load
state.save('data/projects/cafi-survey-2026.json')
state = ProjectState.load('data/projects/cafi-survey-2026.json')
```

---

## Dashboard Features

The Streamlit dashboard provides:

1. **Project Overview** - Progress tracking and artifact browser
2. **Orchestrator** - Workflow planning and next steps
3. **Agent Runner** - Execute individual agents
4. **Artifact Viewer** - Browse generated outputs

Access at: `http://localhost:8501`

---

## Configuration

See `config.py` for settings:

```python
# LLM Settings
LLM_DEFAULT_PROVIDER = 'anthropic'  # or 'openai', 'mock'
CLAUDE_MODEL = 'claude-sonnet-4-5-20250514'
DEFAULT_TEMPERATURE = 0.7
DEFAULT_MAX_TOKENS = 4096

# Directories
PROJECT_DATA_DIR = 'data/projects'
OUTPUTS_DIR = 'outputs'
```

---

## Testing

Test the agent system:

```python
# Test basic import
from agents import get_agent, list_agents
print(list_agents())

# Test agent execution (mock mode)
from config import LLM_DEFAULT_PROVIDER
agent = get_agent('research_prd')
outputs = agent.execute(
    inputs={'project_brief': 'Test', 'data_dictionary': 'Test'},
    context={'use_llm': False}  # Use template fallback
)
print(outputs.keys())
```

---

## Benefits for CAFI Project

1. **Structured Workflow** - Clear progression from PRD to submission
2. **Consistency** - Standardized outputs across analyses
3. **Documentation** - Automatic generation of methods, results
4. **Reproducibility** - State tracking for all artifacts
5. **Synthesis** - AI-assisted interpretation of R outputs

---

## File Structure

```
CAFI-Survey-2026/
├── agents/
│   ├── __init__.py
│   ├── all_agents.py
│   ├── orchestrator_agent.py
│   └── research_prd_agent.py
├── core/
│   ├── __init__.py
│   ├── base_agent.py
│   ├── llm_interface.py
│   └── project_state.py
├── dashboard/
│   ├── __init__.py
│   ├── app.py
│   ├── improved_app.py
│   └── settings_page.py
├── data/
│   ├── projects/          # Agent state files
│   └── *.csv              # Survey data
├── config.py
├── requirements.txt
└── run.py
```

---

*Agent system ready for CAFI Survey Analysis workflow automation.*

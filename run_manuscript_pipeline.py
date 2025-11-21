#!/usr/bin/env python3
"""
CAFI Survey Manuscript Pipeline

Runs all research workflow agents to generate a complete manuscript
from the analysis outputs.
"""

import os
import sys
from pathlib import Path
from datetime import datetime

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

# Load environment variables
from dotenv import load_dotenv
load_dotenv(project_root / ".env")

from agents import get_agent, list_agents
from agents.r_bridge import RBridge
from agents.cafi_context import (
    CAFI_CONTEXT,
    LITERATURE_CONTEXT,
    MODELING_CONTEXT,
    FIGURE_CONTEXT,
    DATA_QA_CONTEXT,
    EDA_CONTEXT
)
from core.project_state import ProjectState, ProjectStage


def create_output_dir():
    """Create directory for agent outputs."""
    output_dir = project_root / "output" / "manuscript"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def save_output(output_dir: Path, filename: str, content: str):
    """Save agent output to file."""
    path = output_dir / filename
    with open(path, 'w') as f:
        f.write(content)
    print(f"  Saved: {path}")


def run_pipeline():
    """Run the complete manuscript generation pipeline."""
    print("=" * 60)
    print("CAFI Survey Manuscript Pipeline")
    print("=" * 60)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Initialize
    output_dir = create_output_dir()
    bridge = RBridge(project_root)

    # Get analysis results
    print("Loading analysis results...")
    results_summary = bridge.get_all_results_summary()
    methods_data = bridge.get_methods_data()
    results_data = bridge.get_results_data()
    print(f"  Loaded {len(bridge.list_tables())} tables, {len(bridge.list_figures())} figures")
    print()

    # Initialize project state
    state = ProjectState(
        project_id="cafi-survey-2026",
        project_title="CAFI Survey Analysis - Coral-Associated Fauna in Mo'orea"
    )

    # Project brief for all agents
    project_brief = f"""
# CAFI Survey Analysis - Research Brief

## Project Title
Coral-associated fauna communities in Pocillopora corals: patterns, predictors, and network structure in Mo'orea, French Polynesia

## Research Questions
1. What factors determine CAFI community composition and diversity?
2. How do coral size, morphology, and condition affect associated fauna?
3. What is the spatial structure of CAFI distributions?
4. What is the network structure of CAFI co-occurrence patterns?

## Study Design
- Observational survey of 204 Pocillopora coral colonies
- 3 reef sites in Mo'orea: HAU, MAT, MRB
- Summer 2019 (June-August)
- 2,847 CAFI individuals identified to 243 OTUs

## Key Hypotheses
H1: Coral volume positively predicts CAFI abundance (habitat size effect)
H2: Wide-branching corals support higher CAFI diversity than tight-branching
H3: CAFI communities show significant spatial autocorrelation
H4: Network modularity reflects functional groupings of CAFI

## Analysis Completed
- 17 R scripts completed
- 114 figures, 76 tables, 11 RDS objects generated
- GLMMs, PERMANOVA, network analysis, machine learning all done

{CAFI_CONTEXT}
"""

    # Context for all agents
    context = {
        'project_root': str(project_root),
        'use_llm': True,
        'cafi_context': CAFI_CONTEXT
    }

    # ========================================
    # Stage 1: Research PRD Agent
    # ========================================
    print("Stage 1: Research PRD Agent")
    print("-" * 40)

    prd_agent = get_agent('research_prd')
    prd_outputs = prd_agent.execute(
        inputs={
            'project_brief': project_brief,
            'data_dictionary': bridge.get_data_dictionary()
        },
        context=context
    )

    save_output(output_dir, "01_research_prd.md", prd_outputs.get('research_prd', ''))
    state.add_agent_output('research_prd', prd_outputs)
    print()

    # ========================================
    # Stage 2: Deep Literature Agent
    # ========================================
    print("Stage 2: Deep Literature Agent")
    print("-" * 40)

    lit_agent = get_agent('literature')
    lit_inputs = {
        'project_brief': project_brief,
        'research_prd': prd_outputs.get('research_prd', ''),
        'seed_papers': LITERATURE_CONTEXT
    }

    lit_outputs = lit_agent.execute(lit_inputs, context)
    save_output(output_dir, "02_literature_synthesis.md", lit_outputs.get('thematic_synthesis', ''))
    save_output(output_dir, "02_gap_analysis.md", lit_outputs.get('gap_analysis', ''))
    state.add_agent_output('literature', lit_outputs)
    print()

    # ========================================
    # Stage 3: Conceptual Framework Agent
    # ========================================
    print("Stage 3: Conceptual Framework Agent")
    print("-" * 40)

    framework_agent = get_agent('framework')
    framework_outputs = framework_agent.execute(
        inputs={
            'research_prd': prd_outputs.get('research_prd', ''),
            'literature_synthesis': lit_outputs.get('thematic_synthesis', '')
        },
        context=context
    )

    save_output(output_dir, "03_conceptual_framework.md", framework_outputs.get('conceptual_frameworks', ''))
    save_output(output_dir, "03_figure1_design.md", framework_outputs.get('figure1_design', ''))
    state.add_agent_output('framework', framework_outputs)
    print()

    # ========================================
    # Stage 4: Data QA Agent
    # ========================================
    print("Stage 4: Data QA Agent")
    print("-" * 40)

    qa_agent = get_agent('data_qa')
    qa_context = {**context, 'domain_context': DATA_QA_CONTEXT}
    qa_outputs = qa_agent.execute(
        inputs={
            'research_prd': prd_outputs.get('research_prd', ''),
            'data_dictionary': bridge.get_data_dictionary(),
            'data_sample': f"Summary stats: {bridge.get_summary_statistics()}"
        },
        context=qa_context
    )

    save_output(output_dir, "04_qa_playbook.md", qa_outputs.get('qa_playbook', ''))
    state.add_agent_output('data_qa', qa_outputs)
    print()

    # ========================================
    # Stage 5: EDA Agent
    # ========================================
    print("Stage 5: EDA Agent")
    print("-" * 40)

    eda_agent = get_agent('eda')
    eda_context = {**context, 'domain_context': EDA_CONTEXT}
    eda_outputs = eda_agent.execute(
        inputs={
            'research_prd': prd_outputs.get('research_prd', ''),
            'cleaned_data': results_summary,
            'qa_report': qa_outputs.get('qa_playbook', '')
        },
        context=eda_context
    )

    save_output(output_dir, "05_eda_plan.md", eda_outputs.get('eda_plan', ''))
    save_output(output_dir, "05_eda_interpretations.md", eda_outputs.get('interpretations', ''))
    state.add_agent_output('eda', eda_outputs)
    print()

    # ========================================
    # Stage 6: Modeling Agent
    # ========================================
    print("Stage 6: Modeling Agent")
    print("-" * 40)

    modeling_agent = get_agent('modeling')
    modeling_context = {**context, 'domain_context': MODELING_CONTEXT}
    modeling_outputs = modeling_agent.execute(
        inputs={
            'research_prd': prd_outputs.get('research_prd', ''),
            'eda_outputs': eda_outputs.get('eda_plan', ''),
            'journal_expectations': 'Target: Ecology or Coral Reefs journal'
        },
        context=modeling_context
    )

    save_output(output_dir, "06_modeling_plan.md", modeling_outputs.get('modeling_prd', ''))
    save_output(output_dir, "06_model_specifications.md", modeling_outputs.get('model_specifications', ''))
    state.add_agent_output('modeling', modeling_outputs)
    print()

    # ========================================
    # Stage 7: Figure Factory Agent
    # ========================================
    print("Stage 7: Figure Factory Agent")
    print("-" * 40)

    figure_agent = get_agent('figures')
    figure_context = {**context, 'domain_context': FIGURE_CONTEXT}
    figure_outputs = figure_agent.execute(
        inputs={
            'research_prd': prd_outputs.get('research_prd', ''),
            'modeling_plan': modeling_outputs.get('modeling_prd', ''),
            'eda_insights': eda_outputs.get('interpretations', '')
        },
        context=figure_context
    )

    save_output(output_dir, "07_figure_inventory.md", figure_outputs.get('figure_inventory', ''))
    save_output(output_dir, "07_figure_captions.md", figure_outputs.get('caption_drafts', ''))
    state.add_agent_output('figures', figure_outputs)
    print()

    # ========================================
    # Stage 8: Scientific Writer Agent
    # ========================================
    print("Stage 8: Scientific Writer Agent - All Sections")
    print("-" * 40)

    writer_agent = get_agent('writer')

    # Write Introduction
    print("  Writing Introduction...")
    intro_outputs = writer_agent.execute(
        inputs={
            'research_prd': prd_outputs.get('research_prd', ''),
            'literature_synthesis': lit_outputs.get('thematic_synthesis', ''),
            'modeling_plan': modeling_outputs.get('modeling_prd', ''),
            'figure_plan': figure_outputs.get('figure_inventory', '')
        },
        context={**context, 'requested_section': 'introduction'}
    )
    save_output(output_dir, "08_introduction.md", intro_outputs.get('draft_introduction', ''))

    # Write Methods
    print("  Writing Methods...")
    methods_outputs = writer_agent.execute(
        inputs={
            'research_prd': prd_outputs.get('research_prd', ''),
            'literature_synthesis': lit_outputs.get('thematic_synthesis', ''),
            'modeling_plan': modeling_outputs.get('modeling_prd', ''),
            'figure_plan': figure_outputs.get('figure_inventory', '')
        },
        context={**context, 'requested_section': 'methods'}
    )
    save_output(output_dir, "08_methods.md", methods_outputs.get('draft_methods', ''))

    # Write Results
    print("  Writing Results...")
    # Include actual results data
    results_context = {
        **context,
        'requested_section': 'results',
        'actual_results': str(results_data)
    }
    results_outputs = writer_agent.execute(
        inputs={
            'research_prd': prd_outputs.get('research_prd', ''),
            'literature_synthesis': lit_outputs.get('thematic_synthesis', ''),
            'modeling_plan': modeling_outputs.get('modeling_prd', ''),
            'figure_plan': figure_outputs.get('figure_inventory', '')
        },
        context=results_context
    )
    save_output(output_dir, "08_results.md", results_outputs.get('draft_results', ''))

    # Write Discussion
    print("  Writing Discussion...")
    discussion_outputs = writer_agent.execute(
        inputs={
            'research_prd': prd_outputs.get('research_prd', ''),
            'literature_synthesis': lit_outputs.get('thematic_synthesis', ''),
            'modeling_plan': modeling_outputs.get('modeling_prd', ''),
            'figure_plan': figure_outputs.get('figure_inventory', '')
        },
        context={**context, 'requested_section': 'discussion'}
    )
    save_output(output_dir, "08_discussion.md", discussion_outputs.get('draft_discussion', ''))

    state.add_agent_output('writer', {
        'introduction': intro_outputs,
        'methods': methods_outputs,
        'results': results_outputs,
        'discussion': discussion_outputs
    })
    print()

    # ========================================
    # Stage 9: Compile Final Manuscript
    # ========================================
    print("Stage 9: Compiling Final Manuscript")
    print("-" * 40)

    manuscript = f"""# Coral-associated fauna communities in Pocillopora corals: patterns, predictors, and network structure in Mo'orea, French Polynesia

## Abstract

[To be written based on completed sections]

---

{intro_outputs.get('draft_introduction', '# Introduction\n\n[Not generated]')}

---

{methods_outputs.get('draft_methods', '# Methods\n\n[Not generated]')}

---

{results_outputs.get('draft_results', '# Results\n\n[Not generated]')}

---

{discussion_outputs.get('draft_discussion', '# Discussion\n\n[Not generated]')}

---

## Acknowledgments

We thank the Mo'orea LTER program and local collaborators for field support.

## References

[To be compiled from citations]

---

*Generated by CAFI Survey Manuscript Pipeline on {datetime.now().strftime('%Y-%m-%d')}*
"""

    save_output(output_dir, "MANUSCRIPT_DRAFT.md", manuscript)

    # Save project state
    state_path = project_root / "data" / "projects" / "cafi-manuscript-state.json"
    state.save(str(state_path))
    print(f"  Project state saved: {state_path}")
    print()

    # ========================================
    # Summary
    # ========================================
    print("=" * 60)
    print("Pipeline Complete!")
    print("=" * 60)
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    print(f"Output directory: {output_dir}")
    print()
    print("Generated files:")
    for f in sorted(output_dir.glob("*.md")):
        size = f.stat().st_size
        print(f"  - {f.name} ({size:,} bytes)")
    print()
    print("Next steps:")
    print("  1. Review generated sections in output/manuscript/")
    print("  2. Edit and refine content")
    print("  3. Add actual statistical results to Results section")
    print("  4. Compile references")
    print("  5. Format for target journal")


if __name__ == "__main__":
    run_pipeline()

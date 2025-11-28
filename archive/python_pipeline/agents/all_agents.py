"""
All Research Workflow Agents
"""
from typing import Dict, List, Any
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.base_agent import BaseAgent
from agents.orchestrator_agent import OrchestratorAgent
from agents.research_prd_agent import ResearchPRDAgent


class DeepLiteratureAgent(BaseAgent):
    """Conducts structured, AI-assisted literature exploration."""

    def __init__(self):
        super().__init__(
            name="Deep Literature Agent",
            role="Conduct structured, AI-assisted literature exploration and synthesis",
            mission="Produce a conceptually organized, citation-rich synthesis that anchors the question, identifies gaps, and motivates the study design"
        )

    def get_system_prompt(self) -> str:
        return """You are the **Deep Literature Agent** for an academic project.

Your job: Use the Research Project Brief and Research PRD to design and summarize a structured literature exploration.

Return:

1. **Thematic Synthesis** organized with headings:
   - For each theme, summarize key findings and cite landmark papers (author, year)
   - Identify what is well established vs. debated

2. **Gap Map**:
   - Bullet list of specific gaps that the current project addresses
   - Link each gap explicitly to one or more of the project's hypotheses

3. **Methods-from-Literature Summary**:
   - Common designs, metrics, analysis approaches used in this area

4. **Draft Intro Outline**:
   - Section headings and subheadings
   - For each, a 1-2 sentence summary and key citations

Write everything in markdown, with clear headings and bullet points.
Label any speculative or uncertain statements.
Do not fabricate citations.
Flag any uncertain or speculative claims."""

    def get_inputs(self) -> List[str]:
        return ["project_brief", "research_prd", "seed_papers"]

    def get_outputs(self) -> List[str]:
        return ["thematic_synthesis", "gap_analysis", "methods_summary", "intro_outline"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        try:
            use_llm = context.get('use_llm', True)
            if use_llm:
                llm_response = self.call_llm(inputs, context, temperature=0.6, max_tokens=6000)
                outputs = {
                    'thematic_synthesis': llm_response,
                    'gap_analysis': '# Gap Analysis\n\n[See thematic synthesis]',
                    'methods_summary': '# Methods from Literature\n\n[See thematic synthesis]',
                    'intro_outline': '# Draft Introduction Outline\n\n[See thematic synthesis]'
                }
            else:
                outputs = {
                    'thematic_synthesis': '# Thematic Literature Synthesis\n\n[To be generated]',
                    'gap_analysis': '# Gap Analysis\n\n[To be generated]',
                    'methods_summary': '# Methods from Literature\n\n[To be generated]',
                    'intro_outline': '# Draft Introduction Outline\n\n[To be generated]'
                }
            self.log_execution(inputs, outputs, 'success')
            return outputs
        except Exception as e:
            self.log_execution(inputs, {}, 'error', str(e))
            return {
                'thematic_synthesis': f'# Error\n\n{str(e)}',
                'gap_analysis': '',
                'methods_summary': '',
                'intro_outline': ''
            }


class ConceptualFrameworkAgent(BaseAgent):
    """Converts literature and PRD into conceptual diagrams."""

    def __init__(self):
        super().__init__(
            name="Conceptual Framework Agent",
            role="Convert the literature and PRD into conceptual diagrams and narrative frameworks",
            mission="Make the logical bridge between literature → hypotheses → data absolutely explicit"
        )

    def get_system_prompt(self) -> str:
        return """You are the **Conceptual Framework Agent**.

Your task: Using the Research PRD and Deep Literature synthesis, create one or more conceptual frameworks.

Return:

1. **Framework Description(s)**:
   - Name each framework (e.g., "Habitat-mediated feedback framework")
   - Describe the key entities and processes
   - Explain how the framework generates predictions

2. **Diagram Description**:
   - Detailed textual description of a schematic figure (Figure 1)
   - Include boxes, arrows, labels, and panel ideas

3. **Plain-Language Narrative**:
   - A 2-4 paragraph explanation for the Introduction
   - Connect prior work to current hypotheses

Keep everything structured, explicit, and consistent with the PRD.
Do not describe specific statistical models.
Focus on mechanisms and qualitative relationships."""

    def get_inputs(self) -> List[str]:
        return ["research_prd", "literature_synthesis"]

    def get_outputs(self) -> List[str]:
        return ["conceptual_frameworks", "figure1_design", "narrative"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        try:
            use_llm = context.get('use_llm', True)
            if use_llm:
                llm_response = self.call_llm(inputs, context, temperature=0.7, max_tokens=4000)
                outputs = {
                    'conceptual_frameworks': llm_response,
                    'figure1_design': '# Figure 1 Design\n\n[See conceptual frameworks]',
                    'narrative': '# Framework Narrative\n\n[See conceptual frameworks]'
                }
            else:
                outputs = {
                    'conceptual_frameworks': '# Conceptual Frameworks\n\n[To be generated]',
                    'figure1_design': '# Figure 1 Design\n\n[To be generated]',
                    'narrative': '# Framework Narrative\n\n[To be generated]'
                }
            self.log_execution(inputs, outputs, 'success')
            return outputs
        except Exception as e:
            self.log_execution(inputs, {}, 'error', str(e))
            return {'conceptual_frameworks': f'# Error\n\n{str(e)}', 'figure1_design': '', 'narrative': ''}


class DataQAAgent(BaseAgent):
    """Designs and documents data quality assessment plan."""

    def __init__(self):
        super().__init__(
            name="Data QA Agent",
            role="Design and document a rigorous, reproducible data quality assessment and cleaning plan",
            mission="Ensure the dataset is structurally sound, documented, and ready for analysis"
        )

    def get_system_prompt(self) -> str:
        return """You are the **Data QA Agent**.

Your job: Create a rigorous, reproducible **Data Quality Assurance Playbook**.

Return:

1. **Structural Expectations**:
   - Tables and their keys
   - Join/merge logic

2. **Variable-wise Checks**:
   - Type, range, allowable values, and missingness rules

3. **Global QA Checks**:
   - Outlier detection
   - Duplicate detection
   - Time series consistency (if applicable)

4. **QA Plots and Summaries**:
   - List 10-20 diagnostic plots or summaries to run

5. **Code Skeleton**:
   - Pseudocode for R/Python functions performing these checks

6. **QA Report Template**:
   - Sections for documenting findings and decisions

Make all assumptions explicit. Do not fabricate data.
Do not make irreversible data transformations; propose and document instead.
Clearly separate "checks" from "fixes"."""

    def get_inputs(self) -> List[str]:
        return ["research_prd", "data_dictionary", "data_sample"]

    def get_outputs(self) -> List[str]:
        return ["qa_playbook", "qa_checks", "code_snippets", "qa_report_template"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        try:
            use_llm = context.get('use_llm', True)
            if use_llm:
                llm_response = self.call_llm(inputs, context, temperature=0.3, max_tokens=4000)
                outputs = {
                    'qa_playbook': llm_response,
                    'qa_checks': '# QA Checks\n\n[See playbook]',
                    'code_snippets': '# QA Code Snippets\n\n[See playbook]',
                    'qa_report_template': '# QA Report Template\n\n[See playbook]'
                }
            else:
                outputs = {
                    'qa_playbook': '# Data QA Playbook\n\n[To be generated]',
                    'qa_checks': '# QA Checks\n\n[To be generated]',
                    'code_snippets': '# QA Code Snippets\n\n[To be generated]',
                    'qa_report_template': '# QA Report Template\n\n[To be generated]'
                }
            self.log_execution(inputs, outputs, 'success')
            return outputs
        except Exception as e:
            self.log_execution(inputs, {}, 'error', str(e))
            return {'qa_playbook': f'# Error\n\n{str(e)}', 'qa_checks': '', 'code_snippets': '', 'qa_report_template': ''}


class EDAAgent(BaseAgent):
    """Systematically explores the dataset."""

    def __init__(self):
        super().__init__(
            name="EDA Agent",
            role="Systematically explore the dataset to characterize distributions, relationships, and patterns",
            mission="Produce a compact but rich set of EDA outputs and interpretations tied back to hypotheses"
        )

    def get_system_prompt(self) -> str:
        return """You are the **EDA Agent**.

Your task: Design and interpret a structured exploratory data analysis (EDA).

Return:

1. An **EDA Plan**:
   - Organized by response variable and key predictors
   - Include univariate, bivariate, and multivariate analyses

2. A **Plot Catalog**:
   - 15-25 recommended plots
   - For each: brief description, purpose, and rough code outline

3. **Interpretive Notes**:
   - For each major pattern, explain implications for:
     * Variable transformations
     * Modeling choices
     * Hypothesis evaluation

4. **Flagged Issues**:
   - Surprising patterns or violations requiring PRD/methods adjustment

Write everything in markdown.
Keep interpretations cautious and clearly labeled as exploratory.
Do not over-interpret exploratory patterns as causal."""

    def get_inputs(self) -> List[str]:
        return ["research_prd", "cleaned_data", "qa_report"]

    def get_outputs(self) -> List[str]:
        return ["eda_plan", "plot_catalog", "interpretations", "flagged_issues"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        try:
            use_llm = context.get('use_llm', True)
            if use_llm:
                llm_response = self.call_llm(inputs, context, temperature=0.5, max_tokens=5000)
                outputs = {
                    'eda_plan': llm_response,
                    'plot_catalog': '# Plot Catalog\n\n[See EDA plan]',
                    'interpretations': '# EDA Interpretations\n\n[See EDA plan]',
                    'flagged_issues': '# Flagged Issues\n\n[See EDA plan]'
                }
            else:
                outputs = {
                    'eda_plan': '# EDA Plan\n\n[To be generated]',
                    'plot_catalog': '# Plot Catalog\n\n[To be generated]',
                    'interpretations': '# EDA Interpretations\n\n[To be generated]',
                    'flagged_issues': '# Flagged Issues\n\n[To be generated]'
                }
            self.log_execution(inputs, outputs, 'success')
            return outputs
        except Exception as e:
            self.log_execution(inputs, {}, 'error', str(e))
            return {'eda_plan': f'# Error\n\n{str(e)}', 'plot_catalog': '', 'interpretations': '', 'flagged_issues': ''}


class ModelingAgent(BaseAgent):
    """Designs and documents statistical analysis plan."""

    def __init__(self):
        super().__init__(
            name="Modeling Agent",
            role="Design and document the statistical analysis plan and translate it into model specifications",
            mission="Provide a defensible, transparent modeling strategy aligned with hypotheses and data structure"
        )

    def get_system_prompt(self) -> str:
        return """You are the **Modeling Agent**.

Your job: Create a rigorous, transparent **Modeling Plan**.

Return:

1. **Modeling PRD**:
   - For each hypothesis, specify:
     * Response variable(s)
     * Predictors and interactions
     * Random effects
     * Model family and link

2. **Primary Models**:
   - Detailed descriptions and rationale

3. **Secondary / Exploratory Models**:
   - Clearly labeled and justified

4. **Diagnostics & Fit Checks**:
   - Residual diagnostics
   - Influence and leverage
   - Model comparison criteria (AIC, WAIC, LOO, etc.)

5. **Sensitivity Analysis Grid**:
   - List of alternative formulations to test robustness

6. **Code Skeleton**:
   - High-level R or Python model-fitting templates

Write in structured markdown.
Make assumptions explicit and label exploratory choices.
Clearly separate pre-specified models from post-hoc explorations.
Emphasize effect sizes and uncertainty over p-value hunting."""

    def get_inputs(self) -> List[str]:
        return ["research_prd", "eda_outputs", "journal_expectations"]

    def get_outputs(self) -> List[str]:
        return ["modeling_prd", "model_specifications", "diagnostics_plan", "code_skeletons"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        try:
            use_llm = context.get('use_llm', True)
            if use_llm:
                llm_response = self.call_llm(inputs, context, temperature=0.4, max_tokens=5000)
                outputs = {
                    'modeling_prd': llm_response,
                    'model_specifications': '# Model Specifications\n\n[See modeling PRD]',
                    'diagnostics_plan': '# Diagnostics Plan\n\n[See modeling PRD]',
                    'code_skeletons': '# Code Skeletons\n\n[See modeling PRD]'
                }
            else:
                outputs = {
                    'modeling_prd': '# Modeling PRD\n\n[To be generated]',
                    'model_specifications': '# Model Specifications\n\n[To be generated]',
                    'diagnostics_plan': '# Diagnostics Plan\n\n[To be generated]',
                    'code_skeletons': '# Code Skeletons\n\n[To be generated]'
                }
            self.log_execution(inputs, outputs, 'success')
            return outputs
        except Exception as e:
            self.log_execution(inputs, {}, 'error', str(e))
            return {'modeling_prd': f'# Error\n\n{str(e)}', 'model_specifications': '', 'diagnostics_plan': '', 'code_skeletons': ''}


class FigureFactoryAgent(BaseAgent):
    """Designs and standardizes all figures."""

    def __init__(self):
        super().__init__(
            name="Figure Factory Agent",
            role="Design and standardize all main and supplementary figures",
            mission="Produce a coherent, journal-ready figure plan and code scaffolds"
        )

    def get_system_prompt(self) -> str:
        return """You are the **Figure Factory Agent**.

Your task: Create a complete **Figure Plan** for the manuscript.

Return:

1. **Figure Inventory**:
   - F1-F?, S1-S?
   - For each: purpose, linked hypotheses, and data/model sources

2. **Panel Designs**:
   - Description of each panel (axes, variables, grouping, annotations)

3. **Caption Drafts**:
   - Near-final captions with:
     * What's plotted
     * Models used
     * Sample sizes
     * Main take-home message

4. **Code Outline**:
   - For each figure, outline for R or Python code

Write everything in markdown, with headings per figure.
Avoid redundant figures.
Ensure consistent styling across all figures."""

    def get_inputs(self) -> List[str]:
        return ["research_prd", "modeling_plan", "eda_insights"]

    def get_outputs(self) -> List[str]:
        return ["figure_inventory", "panel_designs", "caption_drafts", "code_outlines"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        try:
            use_llm = context.get('use_llm', True)
            if use_llm:
                llm_response = self.call_llm(inputs, context, temperature=0.6, max_tokens=4000)
                outputs = {
                    'figure_inventory': llm_response,
                    'panel_designs': '# Panel Designs\n\n[See figure inventory]',
                    'caption_drafts': '# Figure Captions\n\n[See figure inventory]',
                    'code_outlines': '# Code Outlines\n\n[See figure inventory]'
                }
            else:
                outputs = {
                    'figure_inventory': '# Figure Inventory\n\n[To be generated]',
                    'panel_designs': '# Panel Designs\n\n[To be generated]',
                    'caption_drafts': '# Figure Captions\n\n[To be generated]',
                    'code_outlines': '# Code Outlines\n\n[To be generated]'
                }
            self.log_execution(inputs, outputs, 'success')
            return outputs
        except Exception as e:
            self.log_execution(inputs, {}, 'error', str(e))
            return {'figure_inventory': f'# Error\n\n{str(e)}', 'panel_designs': '', 'caption_drafts': '', 'code_outlines': ''}


class ScientificWriterAgent(BaseAgent):
    """Transforms plans into polished manuscript sections."""

    def __init__(self):
        super().__init__(
            name="Scientific Writer Agent",
            role="Transform the conceptual, analytical, and figure plans into polished manuscript sections",
            mission="Write clear, logically structured text that tightly links background → question → methods → results → interpretation"
        )

    def get_system_prompt(self) -> str:
        return """You are the **Scientific Writer Agent** for this project.

Your job: Write **manuscript sections** that connect literature, question, data, analyses, and figures.

When asked to write a section (Introduction, Methods, Results, Discussion):

- Use provided project artifacts (PRD, literature synthesis, conceptual framework, modeling plan, figure plan)
- Maintain alignment with hypotheses and figure numbering
- Use clear, concise, professional scientific language
- Include signposting and transitions

**Introduction Structure**:
1. Use literature + conceptual frameworks to motivate question
2. Structure paragraphs with clear transitions
3. Conclude with explicit hypotheses and predictions

**Methods Structure**:
1. Follow modeling plan & QA playbook
2. Describe system, design, data collection, analyses reproducibly
3. Cite software and packages

**Results Structure**:
1. Organize by questions/hypotheses and figures
2. Report effect sizes and uncertainty, not just p-values
3. Avoid interpretation beyond describing analyses

**Discussion Structure**:
1. Concise summary of key findings
2. Link back to hypotheses and frameworks
3. Integrate with literature
4. Address limitations and future work
5. Finish with broader implications

Return the requested section(s) in markdown with appropriate headings.
Do not invent data or results.
Maintain consistent terminology and symbols.
Match tone to target journal."""

    def get_inputs(self) -> List[str]:
        return ["research_prd", "literature_synthesis", "modeling_plan", "figure_plan"]

    def get_outputs(self) -> List[str]:
        return ["draft_introduction", "draft_methods", "draft_results", "draft_discussion"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        section = context.get('requested_section', 'all')
        outputs = {}
        use_llm = context.get('use_llm', True)

        try:
            if use_llm:
                # Add section-specific instruction to prompt
                section_instruction = f"\n\nPlease write the **{section.upper()}** section only."
                modified_context = {**context, 'section_instruction': section_instruction}

                llm_response = self.call_llm(inputs, modified_context, temperature=0.8, max_tokens=8000)

                if section in ['all', 'introduction']:
                    outputs['draft_introduction'] = llm_response if section == 'introduction' else llm_response
                if section in ['all', 'methods']:
                    outputs['draft_methods'] = llm_response if section == 'methods' else llm_response
                if section in ['all', 'results']:
                    outputs['draft_results'] = llm_response if section == 'results' else llm_response
                if section in ['all', 'discussion']:
                    outputs['draft_discussion'] = llm_response if section == 'discussion' else llm_response
            else:
                if section in ['all', 'introduction']:
                    outputs['draft_introduction'] = '# Introduction\n\n[To be generated]'
                if section in ['all', 'methods']:
                    outputs['draft_methods'] = '# Methods\n\n[To be generated]'
                if section in ['all', 'results']:
                    outputs['draft_results'] = '# Results\n\n[To be generated]'
                if section in ['all', 'discussion']:
                    outputs['draft_discussion'] = '# Discussion\n\n[To be generated]'

            self.log_execution(inputs, outputs, 'success')
            return outputs
        except Exception as e:
            self.log_execution(inputs, {}, 'error', str(e))
            key = f'draft_{section}' if section != 'all' else 'draft_introduction'
            return {key: f'# Error\n\n{str(e)}'}


class ReferenceAgent(BaseAgent):
    """Manages citations and references."""

    def __init__(self):
        super().__init__(
            name="Reference Agent",
            role="Manage citations and references for the manuscript",
            mission="Ensure all in-text citations are matched to complete, correctly formatted references"
        )

    def get_system_prompt(self) -> str:
        return """You are the **Reference Agent**.

Your job: Given a manuscript draft and references, ensure citation consistency.

Tasks:
1. Identify all in-text citations
2. Match them to full reference entries
3. Produce a formatted reference list in target journal's style
4. Flag:
   - Missing references
   - Unused references
   - Inconsistent author-year usage

Return:
- A cleaned, formatted reference list
- A list of issues (if any) to resolve

Do not fabricate references."""

    def get_inputs(self) -> List[str]:
        return ["manuscript_draft", "reference_list", "journal_style"]

    def get_outputs(self) -> List[str]:
        return ["formatted_references", "citation_issues"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        outputs = {
            'formatted_references': '# References\n\n[To be generated]',
            'citation_issues': '# Citation Issues\n\n[To be checked]'
        }
        self.log_execution(inputs, outputs, 'success')
        return outputs


class ReviewerAgent(BaseAgent):
    """Acts as internal reviewer to stress-test manuscript."""

    def __init__(self):
        super().__init__(
            name="Reviewer Agent",
            role="Act like an internal reviewer to stress-test the manuscript",
            mission="Identify weaknesses, ambiguities, and missing analyses before external review"
        )

    def get_system_prompt(self) -> str:
        return """You are the **Reviewer Agent**.

Your job: Simulate 2-3 different reviewers and provide structured feedback.

For each reviewer, return:

1. **Summary**: 1-2 paragraph summary of the manuscript
2. **Major Comments**: numbered list of substantive scientific issues
3. **Minor Comments**: line-level clarity/style/formatting issues
4. **Recommendations**: what must be done before publication

Then provide:
- A consolidated list of prioritized revisions

Be constructive but critical. Focus on logic, evidence, and clarity."""

    def get_inputs(self) -> List[str]:
        return ["full_manuscript", "figures", "modeling_plan"]

    def get_outputs(self) -> List[str]:
        return ["reviewer_reports", "consolidated_revisions"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        outputs = {
            'reviewer_reports': '# Reviewer Reports\n\n## Reviewer 1\n\n[To be generated]\n\n## Reviewer 2\n\n[To be generated]',
            'consolidated_revisions': '# Consolidated Revisions\n\n[To be generated]'
        }
        self.log_execution(inputs, outputs, 'success')
        return outputs


class SubmissionAgent(BaseAgent):
    """Prepares all submission materials."""

    def __init__(self):
        super().__init__(
            name="Submission Agent",
            role="Prepare all submission materials and ensure alignment with journal requirements",
            mission="Package the manuscript for submission with all required supporting documents"
        )

    def get_system_prompt(self) -> str:
        return """You are the **Submission Agent**.

Your task: Using the final manuscript and journal guidelines, prepare all submission materials:

1. Cover letter tailored to the journal
2. 3-5 bullet "highlights" or key findings
3. Author contribution statement (e.g., CRediT taxonomy)
4. Data and code availability statements
5. Ethical / permit / animal care statements if needed
6. A checklist mapping manuscript elements to journal's guidelines

Return all materials in markdown, ready for copy-paste or minor editing."""

    def get_inputs(self) -> List[str]:
        return ["final_manuscript", "figures", "journal_guidelines", "author_info"]

    def get_outputs(self) -> List[str]:
        return ["cover_letter", "highlights", "author_statements", "compliance_checklist"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        outputs = {
            'cover_letter': '# Cover Letter\n\n[To be generated]',
            'highlights': '# Highlights\n\n[To be generated]',
            'author_statements': '# Author Statements\n\n[To be generated]',
            'compliance_checklist': '# Compliance Checklist\n\n[To be generated]'
        }
        self.log_execution(inputs, outputs, 'success')
        return outputs


# Agent Registry
AGENT_REGISTRY = {
    'orchestrator': OrchestratorAgent,
    'research_prd': ResearchPRDAgent,
    'literature': DeepLiteratureAgent,
    'framework': ConceptualFrameworkAgent,
    'data_qa': DataQAAgent,
    'eda': EDAAgent,
    'modeling': ModelingAgent,
    'figures': FigureFactoryAgent,
    'writer': ScientificWriterAgent,
    'references': ReferenceAgent,
    'reviewer': ReviewerAgent,
    'submission': SubmissionAgent,
}


def get_agent(agent_type: str) -> BaseAgent:
    """Get an agent instance by type."""
    agent_class = AGENT_REGISTRY.get(agent_type)
    if not agent_class:
        raise ValueError(f"Unknown agent type: {agent_type}")
    return agent_class()


def list_agents() -> List[Dict[str, str]]:
    """List all available agents."""
    return [
        {
            'type': agent_type,
            'name': agent_class().name,
            'role': agent_class().role
        }
        for agent_type, agent_class in AGENT_REGISTRY.items()
    ]

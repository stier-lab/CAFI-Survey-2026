"""
Research PRD Agent - Academic Problem Definition
"""
from typing import Dict, List, Any
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.base_agent import BaseAgent


class ResearchPRDAgent(BaseAgent):
    """Converts high-level research question into precise, structured PRD."""

    def __init__(self):
        super().__init__(
            name="Research PRD Agent",
            role="Convert a high-level research question + data description into a precise, structured Research Product Requirements Document (PRD)",
            mission="Make the research problem analytically sharp: clearly defined questions, hypotheses, variables, and evaluation criteria"
        )

    def get_system_prompt(self) -> str:
        return """You are the **Research PRD Agent** for an academic project.

Given:
- The Research Project Brief
- Any available data dictionary or notes

Your job is to write a **Research Product Requirements Document (Research PRD)** that makes the problem analytically precise.

The PRD must include:

1. **Umbrella Question and Nested Questions**
   - Restate the high-level research question
   - List of nested sub-questions (Q1, Q2, etc.)

2. **Hypothesis Tree**
   - For each question, list hypotheses (H1a, H1b, etc.)
   - Describe expected effect directions and mechanisms
   - Make predictions testable and falsifiable

3. **Conceptual Causal Diagram (in words)**
   - Key variables and their relationships
   - Confounders, mediators, and moderators
   - Identify assumed causal pathways

4. **Variable Definitions**
   - Response variables (outcomes)
   - Predictor variables (exposures, treatments)
   - Random/grouping factors
   - Experimental treatments
   - For each, specify: type, scale, expected range

5. **Planned Endpoints and Metrics**
   - What outcomes will be used to evaluate hypotheses
   - Any thresholds or decision criteria
   - Primary vs secondary endpoints

6. **Primary vs Secondary Analyses**
   - Which analyses must be run (confirmatory)
   - Which are exploratory or secondary

7. **Risks and Confounders**
   - Potential biases
   - Design limitations
   - Missing data concerns
   - How to address each

Write in clear, structured markdown with headings and bullet points.
Be explicit about assumptions and note them as "Assumption: …" when needed.
Do not write code.
Do not assume data quality; leave that to Data QA Agent.
Explicitly mark any assumptions that are uncertain."""

    def get_inputs(self) -> List[str]:
        return ["project_brief", "data_dictionary"]

    def get_outputs(self) -> List[str]:
        return ["research_prd", "hypothesis_tree", "variable_definitions", "risk_analysis"]

    def execute(self, inputs: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
        """Execute Research PRD generation."""
        try:
            project_brief = inputs.get('project_brief', '')
            data_dictionary = inputs.get('data_dictionary', 'Not provided')

            # Try LLM-powered generation first
            use_llm = context.get('use_llm', True)

            if use_llm:
                try:
                    llm_response = self.call_llm(inputs, context, temperature=0.7, max_tokens=8000)

                    # Parse LLM response into structured outputs
                    outputs = self._parse_llm_response(llm_response, project_brief, data_dictionary)

                    self.log_execution(inputs, outputs, 'success')
                    return outputs

                except Exception as e:
                    # Fall back to template if LLM fails
                    print(f"LLM call failed ({str(e)}), using template")

            # Fallback: return structured template
            outputs = {
                'research_prd': self._generate_prd_template(project_brief, data_dictionary),
                'hypothesis_tree': self._extract_hypothesis_tree(project_brief),
                'variable_definitions': self._extract_variables(project_brief, data_dictionary),
                'risk_analysis': self._identify_risks(project_brief)
            }

            self.log_execution(inputs, outputs, 'success')
            return outputs

        except Exception as e:
            self.log_execution(inputs, {}, 'error', str(e))
            raise

    def _parse_llm_response(self, response: str, project_brief: str, data_dict: str) -> Dict[str, Any]:
        """Parse LLM response into structured outputs."""
        # If LLM is configured and returned real content, use it
        if not response.startswith('[Mock LLM Response]') and not response.startswith('[LLM not configured'):
            # Return the full response as the main PRD
            # In a more sophisticated version, you'd parse sections
            return {
                'research_prd': response,
                'hypothesis_tree': self._extract_hypothesis_tree(project_brief),
                'variable_definitions': self._extract_variables(project_brief, data_dict),
                'risk_analysis': self._identify_risks(project_brief)
            }
        else:
            # Use templates
            return {
                'research_prd': self._generate_prd_template(project_brief, data_dict),
                'hypothesis_tree': self._extract_hypothesis_tree(project_brief),
                'variable_definitions': self._extract_variables(project_brief, data_dict),
                'risk_analysis': self._identify_risks(project_brief)
            }

    def _generate_prd_template(self, brief: str, data_dict: str) -> str:
        """Generate PRD template structure."""
        return f"""# Research Product Requirements Document (PRD)

## 1. Umbrella Question and Nested Questions

### Umbrella Question
[To be extracted from project brief]

### Nested Questions
- Q1: [Specific question 1]
  - Focus: [What aspect]
- Q2: [Specific question 2]
  - Focus: [What aspect]

## 2. Hypothesis Tree

### Q1 Hypotheses
- H1a: [Hypothesis with expected direction]
  - Mechanism: [Why we expect this]
  - Prediction: [What we should observe]
- H1b: [Alternative hypothesis]
  - Mechanism: [Why we expect this]
  - Prediction: [What we should observe]

### Q2 Hypotheses
- H2a: [Hypothesis with expected direction]
  - Mechanism: [Why we expect this]
  - Prediction: [What we should observe]

## 3. Conceptual Causal Diagram

### Key Variables
- **Response Variables**: [List]
- **Primary Predictors**: [List]
- **Confounders**: [List]
- **Mediators**: [List if applicable]
- **Moderators**: [List if applicable]

### Causal Pathways
1. [Predictor] → [Response]: [Description of relationship]
2. [Confounder] → [Predictor] and [Response]: [Description]

## 4. Variable Definitions

### Response Variables
| Variable | Type | Scale | Expected Range | Description |
|----------|------|-------|----------------|-------------|
| [Var1]   | Continuous | Ratio | 0-100 | [Description] |

### Predictor Variables
| Variable | Type | Scale | Expected Range | Description |
|----------|------|-------|----------------|-------------|
| [Var1]   | Categorical | Nominal | A, B, C | [Description] |

### Random Effects / Grouping
| Variable | Type | Levels | Description |
|----------|------|--------|-------------|
| [Group]  | Factor | ~50 | [Description] |

## 5. Planned Endpoints and Metrics

### Primary Endpoints
- [Endpoint 1]: [How measured, what constitutes support for hypothesis]
- [Endpoint 2]: [How measured, what constitutes support for hypothesis]

### Secondary Endpoints
- [Endpoint 3]: [Description]

### Evaluation Criteria
- **Support for H1a**: [Specific criteria, e.g., effect size > X, p < 0.05]
- **Support for H1b**: [Specific criteria]

## 6. Primary vs Secondary Analyses

### Primary (Confirmatory) Analyses
1. [Analysis name]: Tests H1a by [method]
2. [Analysis name]: Tests H1b by [method]

### Secondary (Exploratory) Analyses
1. [Analysis name]: Explores [pattern]
2. [Analysis name]: Sensitivity check for [assumption]

## 7. Risks and Confounders

### Potential Biases
- **[Bias type]**: [Description, how it affects results]
  - Mitigation: [Strategy]

### Design Limitations
- **[Limitation]**: [Description]
  - Impact: [How this limits inference]
  - Mitigation: [If possible]

### Missing Data Concerns
- **[Variable]**: Expected [X]% missing
  - Mechanism: [MCAR, MAR, MNAR?]
  - Handling: [Approach]

### Assumptions to Verify
- Assumption 1: [Description]
- Assumption 2: [Description]

---

*This PRD should be reviewed and refined before proceeding to data analysis.*
"""

    def _extract_hypothesis_tree(self, brief: str) -> Dict[str, List[str]]:
        """Extract hypotheses from brief."""
        return {
            'umbrella_question': 'To be extracted from brief',
            'nested_questions': [],
            'hypotheses': []
        }

    def _extract_variables(self, brief: str, data_dict: str) -> Dict[str, List[Dict[str, str]]]:
        """Extract variable definitions."""
        return {
            'response_variables': [],
            'predictor_variables': [],
            'random_effects': []
        }

    def _identify_risks(self, brief: str) -> List[Dict[str, str]]:
        """Identify potential risks and confounders."""
        return [
            {'risk': 'To be identified from project brief', 'mitigation': 'To be determined'}
        ]


if __name__ == "__main__":
    agent = ResearchPRDAgent()
    print(f"Agent: {agent.name}")
    print(f"Role: {agent.role}")

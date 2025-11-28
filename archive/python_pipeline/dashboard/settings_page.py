"""
Settings Page - API Token Configuration
"""
import streamlit as st
import os
from pathlib import Path


def save_api_key(provider: str, api_key: str):
    """Save API key to .env file."""
    env_path = Path(__file__).parent.parent / ".env"

    # Read existing .env or create new
    env_vars = {}
    if env_path.exists():
        with open(env_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and '=' in line:
                    key, value = line.split('=', 1)
                    env_vars[key.strip()] = value.strip()

    # Update the key
    key_name = f"{provider.upper()}_API_KEY"
    env_vars[key_name] = api_key

    # Write back to .env
    with open(env_path, 'w') as f:
        f.write("# Research Workflow System - API Keys\n")
        f.write("# These keys enable AI-powered agent responses\n\n")
        for key, value in env_vars.items():
            f.write(f"{key}={value}\n")

    # Set in current environment
    os.environ[key_name] = api_key

    return True


def check_api_key(provider: str) -> tuple[bool, str]:
    """Check if API key is configured."""
    key_name = f"{provider.upper()}_API_KEY"
    api_key = os.environ.get(key_name, '')

    if api_key:
        # Mask the key for display
        masked = api_key[:8] + '...' + api_key[-4:] if len(api_key) > 12 else '***'
        return True, masked
    return False, ''


def settings_page():
    """Settings page for configuring API tokens and preferences."""
    st.title("⚙️ Settings")

    st.markdown("""
        Configure your API keys to enable AI-powered agent responses. Your keys are stored
        securely in a local `.env` file and never transmitted anywhere except to the
        respective AI provider's API.
    """)

    st.divider()

    # ============================================================================
    # LLM Provider Configuration
    # ============================================================================

    st.subheader("🤖 AI Provider Configuration")

    st.markdown("""
        The Research Workflow System can use different AI providers for generating
        research content. Choose and configure your preferred provider below.
    """)

    # Tabs for different providers
    tab1, tab2, tab3 = st.tabs(["Anthropic Claude", "OpenAI GPT-4", "Mock (Testing)"])

    # ===== ANTHROPIC TAB =====
    with tab1:
        st.markdown("### Anthropic Claude Configuration")

        st.markdown("""
            **Recommended** - Claude excels at long-form research content, detailed analysis,
            and structured outputs. Perfect for academic research workflows.

            - Model: Claude 3.5 Sonnet
            - Best for: Research PRDs, literature synthesis, manuscript writing
            - Get your API key: [console.anthropic.com](https://console.anthropic.com)
        """)

        # Check current status
        has_key, masked_key = check_api_key('anthropic')

        if has_key:
            st.success(f"✓ Anthropic API key configured: `{masked_key}`")

            col1, col2 = st.columns([1, 3])
            with col1:
                if st.button("Update Key", key="update_anthropic"):
                    st.session_state.show_anthropic_form = True
            with col2:
                if st.button("Remove Key", key="remove_anthropic"):
                    save_api_key('anthropic', '')
                    st.success("API key removed")
                    st.rerun()
        else:
            st.warning("⚠️ No Anthropic API key configured")
            st.session_state.show_anthropic_form = True

        # Show form if needed
        if st.session_state.get('show_anthropic_form', not has_key):
            with st.form("anthropic_form"):
                st.markdown("#### Enter Your Anthropic API Key")

                api_key = st.text_input(
                    "API Key",
                    type="password",
                    placeholder="sk-ant-...",
                    help="Starts with 'sk-ant-' and found in your Anthropic console"
                )

                col1, col2, col3 = st.columns([1, 1, 2])
                with col1:
                    submit = st.form_submit_button("Save Key", type="primary")
                with col2:
                    cancel = st.form_submit_button("Cancel")

                if submit:
                    if api_key:
                        if save_api_key('anthropic', api_key):
                            st.success("✓ API key saved successfully!")
                            st.session_state.show_anthropic_form = False
                            st.rerun()
                    else:
                        st.error("Please enter an API key")

                if cancel:
                    st.session_state.show_anthropic_form = False
                    st.rerun()

    # ===== OPENAI TAB =====
    with tab2:
        st.markdown("### OpenAI GPT-4 Configuration")

        st.markdown("""
            Use OpenAI's GPT-4 for agent responses. Good for general-purpose content
            generation and analysis.

            - Model: GPT-4
            - Best for: General research tasks, data analysis
            - Get your API key: [platform.openai.com](https://platform.openai.com)
        """)

        # Check current status
        has_key, masked_key = check_api_key('openai')

        if has_key:
            st.success(f"✓ OpenAI API key configured: `{masked_key}`")

            col1, col2 = st.columns([1, 3])
            with col1:
                if st.button("Update Key", key="update_openai"):
                    st.session_state.show_openai_form = True
            with col2:
                if st.button("Remove Key", key="remove_openai"):
                    save_api_key('openai', '')
                    st.success("API key removed")
                    st.rerun()
        else:
            st.warning("⚠️ No OpenAI API key configured")
            st.session_state.show_openai_form = True

        # Show form if needed
        if st.session_state.get('show_openai_form', not has_key):
            with st.form("openai_form"):
                st.markdown("#### Enter Your OpenAI API Key")

                api_key = st.text_input(
                    "API Key",
                    type="password",
                    placeholder="sk-...",
                    help="Starts with 'sk-' and found in your OpenAI dashboard"
                )

                col1, col2, col3 = st.columns([1, 1, 2])
                with col1:
                    submit = st.form_submit_button("Save Key", type="primary")
                with col2:
                    cancel = st.form_submit_button("Cancel")

                if submit:
                    if api_key:
                        if save_api_key('openai', api_key):
                            st.success("✓ API key saved successfully!")
                            st.session_state.show_openai_form = False
                            st.rerun()
                    else:
                        st.error("Please enter an API key")

                if cancel:
                    st.session_state.show_openai_form = False
                    st.rerun()

    # ===== MOCK TAB =====
    with tab3:
        st.markdown("### Mock Provider (Testing)")

        st.markdown("""
            The mock provider generates template responses without making API calls.
            Perfect for:

            - Testing the system without API costs
            - Development and debugging
            - Understanding the workflow structure

            **Note:** Mock responses are templates, not AI-generated content.
        """)

        st.info("✓ Mock provider is always available - no configuration needed")

        st.code("""
# Mock provider example output
"[Mock LLM Response]

System: Research PRD Agent prompt...
User: Create Research PRD...

[This is a placeholder response. To get real
AI-generated content, configure an LLM provider.]"
        """, language="text")

    st.divider()

    # ============================================================================
    # Current Configuration Status
    # ============================================================================

    st.subheader("📊 Current Configuration")

    # Detect which provider will be used
    from core.llm_interface import auto_configure_llm

    # Force reconfiguration to pick up new keys
    configured_provider = auto_configure_llm()

    provider_info = {
        'anthropic': {
            'name': 'Anthropic Claude',
            'model': 'Claude 3.5 Sonnet',
            'status': '✓ Active' if configured_provider == 'anthropic' else '○ Available'
        },
        'openai': {
            'name': 'OpenAI GPT-4',
            'model': 'GPT-4',
            'status': '✓ Active' if configured_provider == 'openai' else '○ Available'
        },
        'mock': {
            'name': 'Mock Provider',
            'model': 'Template Responses',
            'status': '✓ Active' if configured_provider == 'mock' else '○ Available'
        }
    }

    current = provider_info[configured_provider]

    st.success(f"""
    **Active Provider:** {current['name']}
    **Model:** {current['model']}
    **Status:** {current['status']}
    """)

    # Show all providers status
    with st.expander("View All Providers"):
        for provider, info in provider_info.items():
            has_key, _ = check_api_key(provider)
            status_icon = "✓" if has_key or provider == 'mock' else "○"
            st.markdown(f"**{status_icon} {info['name']}** - {info['model']}")

    st.divider()

    # ============================================================================
    # Usage Guidelines
    # ============================================================================

    st.subheader("📖 Usage Guidelines")

    with st.expander("API Key Security"):
        st.markdown("""
        ### How Your Keys Are Stored

        - Keys are saved to `.env` file in the project root
        - The `.env` file is in `.gitignore` (never committed to version control)
        - Keys are only used to make API calls to the respective provider
        - Keys are never transmitted to any third party

        ### Best Practices

        1. **Never share your API keys** with anyone
        2. **Don't commit `.env` to git** (already in .gitignore)
        3. **Rotate keys regularly** for security
        4. **Use separate keys** for development and production
        5. **Monitor your API usage** in the provider's dashboard
        """)

    with st.expander("Cost Considerations"):
        st.markdown("""
        ### API Costs

        Both Anthropic and OpenAI charge per token (roughly per word):

        - **Anthropic Claude 3.5 Sonnet:** ~$3 per million input tokens, ~$15 per million output tokens
        - **OpenAI GPT-4:** ~$30 per million input tokens, ~$60 per million output tokens

        ### Typical Research Project Costs

        - Research PRD generation: $0.10 - $0.50
        - Literature synthesis: $0.50 - $2.00
        - Full manuscript draft: $2.00 - $10.00
        - Complete project workflow: $5.00 - $25.00

        ### Cost Optimization Tips

        1. Start with mock provider to understand the workflow
        2. Use Claude (cheaper than GPT-4 with similar quality)
        3. Review and refine prompts before running
        4. Save and reuse good outputs
        5. Use lower temperature (0.3-0.5) for deterministic outputs
        """)

    with st.expander("Troubleshooting"):
        st.markdown("""
        ### Common Issues

        **"API key not found" error**
        - Make sure you saved the key in this settings page
        - Check that `.env` file exists in project root
        - Restart the dashboard after adding keys

        **"Invalid API key" error**
        - Verify the key is correct (copy-paste from provider console)
        - Check for extra spaces before/after the key
        - Ensure the key hasn't been revoked in provider dashboard

        **"Rate limit exceeded" error**
        - You've hit the provider's usage limits
        - Wait a few minutes and try again
        - Upgrade your plan with the provider if needed

        **Mock provider being used instead of configured provider**
        - Restart the dashboard application
        - Check that API key is correctly formatted
        - Try removing and re-adding the key
        """)

    st.divider()

    # ============================================================================
    # Test API Connection
    # ============================================================================

    st.subheader("🧪 Test API Connection")

    st.markdown("Test your configured API key to ensure it's working correctly.")

    if st.button("Test Connection", type="primary"):
        try:
            from core.llm_interface import get_llm

            with st.spinner("Testing API connection..."):
                llm = get_llm()

                response = llm.generate(
                    system_prompt="You are a helpful assistant.",
                    user_prompt="Say 'Hello! API connection successful.' and nothing else.",
                    max_tokens=50,
                    temperature=0.3
                )

                if "API connection successful" in response or not response.startswith("[Mock"):
                    st.success(f"""
                    ✓ API Connection Successful!

                    **Provider:** {llm.provider_name}

                    **Test Response:**
                    {response}
                    """)
                else:
                    st.warning(f"""
                    ⚠️ Using Mock Provider

                    No API key configured. The system is using template responses.
                    Configure an API key above to enable AI-powered responses.
                    """)
        except Exception as e:
            st.error(f"""
            ✗ Connection Failed

            **Error:** {str(e)}

            Please check your API key and try again.
            """)


if __name__ == "__main__":
    settings_page()

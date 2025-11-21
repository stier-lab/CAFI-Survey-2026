#!/usr/bin/env python3
"""
Main entry point for Research Workflow System
"""
import sys
import subprocess
import os


def check_dependencies():
    """Check if required dependencies are installed."""
    try:
        import streamlit
        return True
    except ImportError:
        return False


def install_dependencies():
    """Install required dependencies."""
    print("Installing dependencies...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"])
    print("Dependencies installed successfully!\n")


def run_tests():
    """Run the test suite."""
    print("Running test suite...\n")
    result = subprocess.run([sys.executable, "test_agents.py"])
    return result.returncode == 0


def start_dashboard():
    """Start the Streamlit dashboard."""
    print("\n" + "=" * 60)
    print("STARTING RESEARCH WORKFLOW DASHBOARD")
    print("=" * 60)
    print("\nThe dashboard will open in your browser at: http://localhost:8501")
    print("\nPress Ctrl+C to stop the dashboard\n")

    dashboard_path = os.path.join("dashboard", "app.py")
    subprocess.run([sys.executable, "-m", "streamlit", "run", dashboard_path])


def main():
    """Main entry point."""
    print("=" * 60)
    print("RESEARCH WORKFLOW SYSTEM")
    print("=" * 60)
    print()

    # Check dependencies
    if not check_dependencies():
        print("Required dependencies not found.")
        response = input("Would you like to install them now? (y/n): ")
        if response.lower() == 'y':
            install_dependencies()
        else:
            print("Please run: pip install -r requirements.txt")
            sys.exit(1)

    # Ask if user wants to run tests
    print("\nWould you like to run the test suite first?")
    response = input("(Recommended for first run) (y/n): ")

    if response.lower() == 'y':
        if not run_tests():
            print("\nSome tests failed. You can still proceed, but there may be issues.")
            response = input("Continue anyway? (y/n): ")
            if response.lower() != 'y':
                sys.exit(1)

    # Start dashboard
    start_dashboard()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nDashboard stopped. Goodbye!")
        sys.exit(0)
    except Exception as e:
        print(f"\nError: {str(e)}")
        sys.exit(1)

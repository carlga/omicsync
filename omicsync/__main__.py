"""Entry point for the omicsync package."""

import sys


def main():
    """Execute the main entry point with proper error handling."""
    try:
        from .main import main as omicsync_main
        omicsync_main()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\nFatal error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

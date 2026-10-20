"""Groq API client for the Dyscovery research assistant."""

import os

from dotenv import load_dotenv
from groq import APIConnectionError, APIStatusError, APITimeoutError, Groq


load_dotenv()

DEFAULT_MODEL = "openai/gpt-oss-120b"


class GroqConfigurationError(RuntimeError):
    """Raised when the Groq client is not configured."""


def create_client() -> Groq:
    """Create a Groq client using the local environment configuration."""
    api_key = os.getenv("GROQ_API_KEY")

    if not api_key:
        raise GroqConfigurationError(
            "GROQ_API_KEY is missing. Configure it in your local .env file."
        )

    return Groq(
        api_key=api_key,
        timeout=30.0,
        max_retries=1,
    )


def generate_response(
    prompt: str,
    system_prompt: str = (
        "You are a biological research assistant. "
        "Distinguish established facts from uncertainty. "
        "Do not invent scientific references or citations."
    ),
    model: str | None = None,
) -> str:
    """Generate a response from Groq with diagnostic error handling."""
    if not prompt or not prompt.strip():
        raise ValueError("The prompt must not be empty.")

    client = create_client()
    selected_model = model or os.getenv("GROQ_MODEL", DEFAULT_MODEL)

    try:
        response = client.chat.completions.create(
            model=selected_model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": prompt},
            ],
            temperature=0.2,
            max_tokens=1200,
        )

        if not response.choices:
            raise RuntimeError(
                "Groq returned no completion choices. "
                "Try again or check model availability."
            )

        choice = response.choices[0]
        content = choice.message.content
        if isinstance(content, str) and content.strip():
            return content.strip()

        finish_reason = getattr(choice, "finish_reason", None)
        refusal = getattr(choice.message, "refusal", None)

        if refusal:
            raise RuntimeError(f"Groq declined the request: {refusal}")

        raise RuntimeError(
            "Groq returned an empty response. "
            f"Model={selected_model}; finish_reason={finish_reason!r}. "
            "Try a shorter evidence prompt or retry the search."
        )

    except APITimeoutError as exc:
        raise RuntimeError(
            "Groq request timed out. Check your network and try again."
        ) from exc
    except APIConnectionError as exc:
        raise RuntimeError(
            "Could not connect to Groq. Check your network connection."
        ) from exc
    except APIStatusError as exc:
        raise RuntimeError(
            f"Groq API returned HTTP {exc.status_code}. "
            "Check your API key, model availability, and usage limits."
        ) from exc

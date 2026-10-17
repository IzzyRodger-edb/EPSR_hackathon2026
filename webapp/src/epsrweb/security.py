"""Passwords, session cookies and encryption of stored Overleaf tokens.

Passwords use `hashlib.scrypt` from the standard library, so there is no
bcrypt/argon2 wheel to build. Session cookies are signed, not encrypted: they
carry a user id and an expiry and nothing secret. Overleaf tokens *are* secret
-- one token grants write access to every project of a user's account -- so they
are encrypted at rest with a key derived from `EPSRWEB_SECRET_KEY` and are never
written to a log or rendered into a page.
"""

from __future__ import annotations

import base64
import hashlib
import hmac
import secrets
import time

from cryptography.fernet import Fernet, InvalidToken
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.kdf.hkdf import HKDF

_SCRYPT_N = 2**14
_SCRYPT_R = 8
_SCRYPT_P = 1
_DKLEN = 32

MIN_PASSWORD_LENGTH = 10


# --------------------------------------------------------------------------- passwords


def hash_password(password: str) -> str:
    salt = secrets.token_bytes(16)
    digest = hashlib.scrypt(
        password.encode("utf-8"), salt=salt, n=_SCRYPT_N, r=_SCRYPT_R, p=_SCRYPT_P, dklen=_DKLEN
    )
    return "scrypt${}${}${}${}${}".format(
        _SCRYPT_N,
        _SCRYPT_R,
        _SCRYPT_P,
        base64.b64encode(salt).decode("ascii"),
        base64.b64encode(digest).decode("ascii"),
    )


def verify_password(password: str, stored: str) -> bool:
    try:
        scheme, n, r, p, salt_b64, digest_b64 = stored.split("$")
        if scheme != "scrypt":
            return False
        salt = base64.b64decode(salt_b64)
        expected = base64.b64decode(digest_b64)
        actual = hashlib.scrypt(
            password.encode("utf-8"),
            salt=salt,
            n=int(n),
            r=int(r),
            p=int(p),
            dklen=len(expected),
        )
    except (ValueError, TypeError):
        return False
    return hmac.compare_digest(actual, expected)


def password_problem(password: str) -> str | None:
    """Return a message when the password is unacceptable, else None."""
    if len(password) < MIN_PASSWORD_LENGTH:
        return f"The password must be at least {MIN_PASSWORD_LENGTH} characters."
    if password.strip() == "":
        return "The password cannot be blank."
    return None


# --------------------------------------------------------------- signed session cookies


def _sign(secret: bytes, message: bytes) -> str:
    mac = hmac.new(secret, message, hashlib.sha256).digest()
    return base64.urlsafe_b64encode(mac).decode("ascii").rstrip("=")


def make_session_token(secret: bytes, user_id: int, max_age: int) -> str:
    """`<user id>.<expiry>.<nonce>.<mac>`; the nonce makes tokens distinguishable."""
    expires = int(time.time()) + max_age
    body = f"{user_id}.{expires}.{secrets.token_urlsafe(8)}"
    return f"{body}.{_sign(secret, body.encode('ascii'))}"


def read_session_token(secret: bytes, token: str | None) -> int | None:
    """Return the user id of a valid, unexpired token, else None."""
    if not token:
        return None
    parts = token.rsplit(".", 1)
    if len(parts) != 2:
        return None
    body, mac = parts
    if not hmac.compare_digest(_sign(secret, body.encode("ascii")), mac):
        return None
    fields = body.split(".")
    if len(fields) != 3:
        return None
    try:
        user_id, expires = int(fields[0]), int(fields[1])
    except ValueError:
        return None
    if expires < time.time():
        return None
    return user_id


# ------------------------------------------------------------------------------- CSRF


def make_csrf_token(secret: bytes, session_token: str) -> str:
    """A token bound to the session, so it cannot be replayed into another one."""
    return _sign(secret + b"csrf", session_token.encode("utf-8"))


def check_csrf_token(secret: bytes, session_token: str | None, submitted: str | None) -> bool:
    if not session_token or not submitted:
        return False
    return hmac.compare_digest(make_csrf_token(secret, session_token), submitted)


# ------------------------------------------------------------------ secret-at-rest box


def _fernet(secret: bytes) -> Fernet:
    key = HKDF(
        algorithm=hashes.SHA256(), length=32, salt=None, info=b"epsrweb.overleaf.token.v1"
    ).derive(secret)
    return Fernet(base64.urlsafe_b64encode(key))


def encrypt_secret(secret: bytes, plaintext: str) -> str:
    return _fernet(secret).encrypt(plaintext.encode("utf-8")).decode("ascii")


def decrypt_secret(secret: bytes, ciphertext: str) -> str | None:
    """Decrypt, or None if the key changed or the value is corrupt."""
    try:
        return _fernet(secret).decrypt(ciphertext.encode("ascii")).decode("utf-8")
    except (InvalidToken, ValueError):
        return None

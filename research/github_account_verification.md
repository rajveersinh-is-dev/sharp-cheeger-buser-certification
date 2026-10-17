# GitHub account verification (integrity record)

Checked: **2026-10-18**, before creating any remote.

## What I checked, and what I found

`gh auth status` reported a discrepancy that had to be resolved before acting:

```
$ gh auth status
  X Logged in to github.com account daveaddams91-dev (keyring)
  - Active account: true
  ...
  X Failed to log in to github.com account Raj123-0 (default)
```

This suggested the active account might be `daveaddams91-dev`.

However, the **live API disagrees**:

```
$ gh api user --jq '{login: .login, id: .id, name: .name, created_at: .created_at}'
{"created_at":"2026-08-21T09:30:37Z","id":319411334,"login":"rajveersinh-is-dev","name":"Rajveersinh Pardeshi"}
```

## Resolution

The **live API response is authoritative** and it says the active account is:

* **login:** `rajveersinh-is-dev`
* **id:** `319411334`
* **display name:** `Rajveersinh Pardeshi`
* **created:** 2026-08-21

The `daveaddams91-dev` string in `gh auth status` is evidently a local
authentication-store entry that does not correspond to the identity the token
actually carries. Anyone re-running `gh auth status` alone would reach the wrong
conclusion; the token scope list includes `repo` and `user`, so `gh api user` is
available and is the correct way to establish identity.

Note also that `Raj123-0` (a second, invalid stored account) is **not** active
and was not used.

## Consequence for this task

The repository is created under **`rajveersinh-is-dev`**, matching the
token-holder's actual identity and display name, which also matches the author
name required for this task (**Rajveersinh Vishal Pardeshi**).

The repository name `sharp-cheeger-buser-certification` was checked for
availability under that account and **does not already exist**:

```
$ gh repo view rajveersinh-is-dev/sharp-cheeger-buser-certification
GraphQL: Could not resolve to a Repository with the name ...
```

(The identical check against `daveaddams91-dev` also reported "not found", so the
name is free either way; ownership follows the verified active account.)

## What this does *not* assert

This record asserts only which GitHub account the CLI is authenticated as. It
does **not** assert any real-world affiliation, employment, institutional
relationship, or identity beyond what the GitHub API returned. No institutional
affiliation is claimed anywhere in this project; the author is designated
"Independent Researcher" throughout, which is a statement about how the work is
presented and not a claim of employment.
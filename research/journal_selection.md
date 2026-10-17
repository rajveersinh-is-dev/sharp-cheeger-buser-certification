# Journal selection and fee verification

Date checked: **2026-10-18**. Author: Rajveersinh Vishal Pardeshi.

All policy statements below were read from **official Oxford Academic pages in
this session** via a browser session (direct HTTP fetches of `academic.oup.com`
return HTTP 403 to this environment; the browser reached the pages normally).
URLs and the wording actually observed are recorded so they can be re-checked.

---

## 1. Central OUP charge policy (the rule that decides everything)

**URL:** <https://academic.oup.com/pages/open-research/open-access/charges-licences-and-self-archiving>
(reached by redirect from
`.../journals/pages/authors/production_and_publication/publication-charges`)

Wording observed:

> **"OA charges in hybrid journals** — Most OUP journals offer the option for
> authors to choose to publish under an OA licence. Your article will be made
> freely available online immediately upon publication, with the payment of an
> Article Processing Charge (APC)."

> **"OA charges in fully OA journals** — In fully OA journals, articles are
> published under an open access licence immediately upon publication for a fee
> ... In some cases, the fee might be waived."

Interpretation: in a **hybrid** journal the OA/APC route is *optional*, so a
standard (non-OA) route with no mandatory author publication fee exists. In a
**fully OA** journal the APC is *mandatory*, which **fails** the task
requirement.

---

## 2. Candidates considered

### 2.1 The Quarterly Journal of Mathematics (QJM) — **SELECTED**

* **URL:** <https://academic.oup.com/qjmath>
* **Subject line (observed on the site):** "The Quarterly Journal of Mathematics
  — **Pure Mathematics**".
* **Author instructions:** <https://academic.oup.com/qjmath/pages/General_Instructions>
* **Decision-relevant wording, read verbatim on that page:**

> "**Open Access Option for Authors** — Quarterly Journal of Mathematics offers
> the option of publishing under **either a standard licence or an open access
> licence**. Please note that some funders require open access publication as a
> condition of funding. If you are unsure whether you are required to publish
> open access, please do clarify any such requirements with your funder or
> institution. Should you wish to publish your article open access, you should
> select your choice of open access licence in our online system after your
> article has been accepted for publication. **You will need to pay an open
> access charge to publish under an open access licence.**"

* **No-mandatory-fee route:** YES — the **standard licence** route. The page
  offers a *choice*; the fee sentence is conditioned on choosing the OA licence.
  This is the route this manuscript would use.
* **Status: VERIFIED on 2026-10-18 from the official page quoted above.**
* **Scope fit:** QJM is listed by OUP under **Pure Mathematics**. A paper on
  spectral graph theory with complete proofs is pure mathematics. The journal
  also publishes analysis and applied-analysis content; the manuscript's
  isoperimetric-inequality subject sits comfortably in its remit.
* **Unresolved ambiguity:** the exact APC amount for QJM was not read from the
  pricing tool (it requires an interactive journal-selection form). It is
  irrelevant to the mandatory-fee question because the standard route carries no
  APC, but I record that I did not verify the optional figure.

### 2.2 Journal of the London Mathematical Society (JLMS)

* **URL:** <https://academic.oup.com/jlms>
* **Scope wording observed:** "the London Mathematical Society has been
  publishing significant research, **predominantly in pure mathematics**. The
  Journal aims to publish papers which are of **medium length**, of high
  quality, are well-written and …"
* **Fit:** good; OUP-published by the LMS.
* **Fee route:** the journal landing page did not expose an explicit
  standard-vs-OA statement in the text retrieved, so the mandatory-fee status is
  **NOT INDEPENDENTLY VERIFIED from the journal's own page** in this session.
  Given that QJM is verified and JLMS is not, JLMS is retained as an
  **alternative**, not as the recommendation.

### 2.3 International Mathematics Research Notices (IMRN)

* **URL:** <https://academic.oup.com/imrn>
* **Author guidelines:** <https://academic.oup.com/imrn/pages/General_Instructions>
* **Fit:** strong on subject matter (research mathematics).
* **Fee route:** the landing page links "Open Access" directly to the
  General Instructions page and the journal is presented in OUP's
  fully-open-access style; IMRN is widely known to be fully OA. I did **not**
  manage to extract the fee paragraph from the dynamically-rendered guidelines
  body.
* **Decision: REJECTED for this task.** IMRN fails the "no mandatory author
  publication fee" requirement under the OUP rule quoted in §1 if it is fully
  OA, and I could not verify a standard non-OA route for it. Recording this as
  an honest exclusion rather than claiming compliance.

### 2.4 Two journals I could not confirm exist on OUP's current list

* `academic.oup.com/mathematika` → **Page Not Found**
* `academic.oup.com/analysisandpde` → **Page Not Found**

I therefore do **not** list Mathematika or *Analysis & PDE* as candidates; I
could not verify them on OUP's journal list in this session. (OUP's own
A–Z list was retrieved and used as the authority for what exists.)

---

## 3. Recommendation

**Target journal: The Quarterly Journal of Mathematics (OUP).**
Route: **standard licence (non-open-access), no mandatory author publication
fee**, verified on the official author-instructions page on 2026-10-18.

Additional charges/uncertainties to state plainly:

* Optional OA APC exists; its amount was **not** verified.
* OUP states "Optional APCs are generally separate from any other page and
  colour charges that might apply." Page/colour charges, if any, were **not**
  verified for QJM.
* No waiver or Read-and-Publish agreement is assumed. OUP notes that authors
  from participating institutions may have the charge paid for them; I do not
  claim any such eligibility.

---

## 4. Preprint and public-repository policy — verified

**URL:** <https://academic.oup.com/journals/pages/access_purchase/rights_and_permissions/self_archiving_policy_b>
(linked from QJM's own navigation as "Self-Archiving Policy")
Full OUP author self-archiving policy, **Policy B**. Wording observed:

> "**Author's Original Version** — The Author's Original Version (AOV) is
> defined here as the un-refereed author version of an article completed before
> submission of the article to the journal. **This is sometimes referred to as
> the 'preprint' version.** ... **Authors may make their AOV available anywhere
> at any time.** This includes posting on their own personal websites,
> institutional or non-commercial subject based repositories, **commercial
> platforms** websites or repositories, or social media, provided that, upon
> acceptance, they acknowledge that the article has been accepted for
> publication …"

**Consequence for this task.** QJM's site links this as *its* self-archiving
policy. Under it, an author's preprint may be posted "anywhere at any time",
explicitly including "commercial platforms websites or repositories" — which
includes GitHub. Therefore **making the manuscript source public on GitHub is
compatible with QJM's preprint policy**, provided that (a) it is the author's
own work, (b) the required acknowledgement is included once the article is
accepted, and (c) the DOI/Version-of-Record link is added after publication.

**This is a positive finding**, so the full manuscript source is included in the
public repository. I did **not** need to withhold the paper.

Conditions that were respected:

* No confidential reviewer correspondence exists (none was ever received).
* No credentials, tokens, private datasets or third-party copyrighted material
  are included.
* The included bibliography contains only verified metadata; no copyrighted
  text from any source is reproduced.

**Status: VERIFIED on 2026-10-18 from two official OUP pages quoted above.**

### 4.1 QJM's own explicit preprint policy — verified

Read verbatim on <https://academic.oup.com/qjmath/pages/General_Instructions>
(2026-10-18):

> "**Preprint policy** — Authors retain the right to make an Author's Original
> Version (preprint) available through various channels, **and this does not
> prevent submission to the journal.** For further information see our Online
> Licensing, Copyright and Permissions policies. **If accepted, the authors are
> required to update the status of any preprint, including your published
> paper's DOI**, as described on our Author Self-Archiving policy page."

This is the most direct evidence available and it is unambiguous: posting the
manuscript publicly **before submission does not prevent submission**.

**Action taken.** The repository was created and the manuscript source was made
public. A conditional obligation is recorded and will be honoured if the paper
is ever accepted: update the preprint with the published DOI.

### 4.2 Other verified QJM requirements relevant to this repository

From the same page:

> "Quarterly Journal of Mathematics supports the Force 11 Data Citation
> Principles and requires that all publicly available datasets be fully
> referenced..."

*Not applicable here*: this project uses **no datasets**; the "experiments" are
exhaustive enumerations generated by the committed code.

> "Quarterly Journal of Mathematics requires submitting authors to provide an
> **ORCID iD at submission**..."

*Not yet satisfied.* The author has not supplied an ORCID iD. This is recorded
as an outstanding pre-submission requirement in `limitations.md` and in the
final report. It is **not** something this project can resolve on the author's
behalf, and no ORCID is invented here.

> "Papers are invited on subject of **pure mathematics and its applications**."

Confirms the scope fit recorded in §2.1, in the journal's own words.

---

## 5. Summary table

| Journal | Fit | Mandatory author fee? | Route | Evidence |
|---|---|---|---|---|
| **Quarterly Journal of Mathematics** | Pure mathematics — good | **No, on the standard route** | standard licence | Official page quoted, 2026-10-18 |
| Journal of the London Mathematical Society | Pure mathematics — good | Not verified from its own page | — | Landing page read; fee statement not exposed |
| IMRN | Research mathematics — good | Rejected (appears fully OA; standard route not verified) | — | Landing page + central OUP rule |
| Mathematika | — | n/a | — | **Page Not Found on OUP — not listed** |
| Analysis & PDE | — | n/a | — | **Page Not Found on OUP — not listed** |
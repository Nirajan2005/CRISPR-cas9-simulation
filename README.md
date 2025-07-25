 CRISPR-Cas9 Therapeutic Gene Editing Simulator  
**By Nirajan Raj Kunwar** | 19 y/o Student Researcher, Nepal 🇳🇵  


A Python-based interactive simulator that models CRISPR-Cas9 gene editing across human organs, with delivery vector dynamics, mutation types, immune clearance, and on/off-target editing logic.  
Built for research learning, science outreach, and therapeutic insight.

---

 What This Simulates

This project simulates how gene editing works in a real-world context:

- **Organ-Specific Delivery** – Pick a human organ system for therapy  
- **Vector Choice** – Choose between LNPs, retroviruses, or adenoviruses  
- **Immune Clearance Risk** – Simulates vector success probability  
- **CRISPR Guide RNA Matching** – Detects on-target vs off-target edits  
- **Edit Type Logic** – Insert, Replace, Silence mutations modeled  
- **DNA Repair Simulation** – Supports NHEJ & HDR outcomes

---

 Real-World Relevance

This sim mirrors actual challenges in modern gene therapy:

| Feature | Real-World Counterpart |
|--------|------------------------|
| LNPs | mRNA vaccine & CRISPR delivery (e.g. COVID, ex vivo therapy) |
| Viral vectors | Retro/Adenoviruses used in CAR-T & gene therapy |
| DNA Repair | NHEJ (knockouts) & HDR (precise edits) |
| Off-Target Edits | Key issue in CRISPR safety research |
| Organ Delivery | Major hurdle in in vivo therapies |

---

 Run It Yourself

1. Clone the repo:
   ```bash
   git clone https://github.com/Nirajan2005/CRISPR-cas9-simulation
   cd CRISPR-cas9-simulation
python crispr_delivery_simulation.py

Choose target organ:
1. Liver
2. Muscle
3. Brain
...
Enter your gRNA: AGCTTAGGCTA...
Edit type: Insert
Strand: sense

Why I Built This
I’m Nirajan — a 19-year-old from Nepal. I built this to help make complex biotech concepts like CRISPR accessible, visual, and exciting, especially for students in low-resource countries.

This is part of my larger project, 100 Brilliant Minds Nepal, to spark a science movement here and empower students through biotech, AI, and creativity.

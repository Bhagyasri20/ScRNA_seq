# 🔍 Python Command-Line Spell Checker

A lightweight, interactive command-line **Spell Checker Application** built in Python. The app processes user-inputted text in real time, identifies spelling mistakes, and outputs corrected alternatives directly to the console.

## 🚀 Features

*   **Interactive CLI Loop:** Runs continuously in the terminal until explicitly closed by typing `exit`.
*   **Step-by-Step Logging:** Prints a detailed breakdown of which misspelled words are being modified and what they are changing into.
*   **Tokenized Validation:** Splits incoming phrases into individual word tokens to systematically verify them against a frequency list.
*   **Smart Levenshtein Corrections:** Powered by `pyspellchecker` to find words with the lowest edit distance.

---

## 🛠️ Tech Stack

*   **Language:** Python 3.12+
*   **Primary Library:** `pyspellchecker`
*   **Interface:** Command Line / Terminal

---

## ⚙️ Installation & Setup

Running by following these simple installation steps.

### Prerequisites
[Python](https://python.org) installed on machine.

### 1. Clone the Repository
```bash
git clone https://github.com
cd YOUR_REPO_NAME
```

### 2. Install the Required Library
Install the `pyspellchecker` package using pip:
```bash
pip install pyspellchecker
```

### 3. Run the Application
Execute the script using Python:
```bash
python spell_checker.py
```

## 💡 How It Works

The core program logic is structured inside the `SpellCheckerApp` class:
1.  **Initialization (`__init__`):** Spawns an instance of the `SpellChecker` module.
2.  **Splitting (`correct_text`):** Breaks down the target sentence into separate words using a space separator.
3.  **Correction Pipeline:** Checks if the lowercase variation matches the dictionary validation; if a typo is found, it logs the change and attaches the correct choice to the output string.
4.  **Reconstruction:** Merges the tracked tokens back together into a readable sentence string.

---

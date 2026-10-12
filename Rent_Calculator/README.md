
Roommate Rent & Expense Calculator

A lightweight Python application designed to calculate and split shared living expenses equally among roommates, flatmates, or hostel residents.

📌 Overview

Managing shared household costs can be tedious. This script automates the process by taking monthly costs—including base rent, collective food orders, and variable electricity consumption—and calculating the exact total amount each person needs to contribute.

✨ Features

• Interactive CLI Inputs: Dynamic prompts for rent, food bills, and power metrics.
• Metered Electricity Breakdown: Computes the total power bill based on exact units consumed and per-unit utility rates.
• Equal Split Evaluation: Aggregates all totals and divides them cleanly by the number of residents using integer floor division.

🛠️ Technical Breakdown


Core Logic

The application processes inputs sequentially using the following calculation steps:
python
# 1. Calculate the total power bill
Total_electricity_bill = Electricity_spend * Charge_per_unit

# 2. Sum all expenses and divide by total residents
Expense = (Food + Rent + Total_electricity_bill) // Persons
Use code with caution.

Input Parameters Required

• Rent: The base cost of the hostel or flat.
• Food: Total cost of collective food ordered for snacking.
• Electricity_spend: Absolute power units consumed during the billing cycle.
• Charge_per_unit: The utility cost rate per single unit of electricity.
• Persons: The total count of individuals living together.

🚀 Getting Started


Prerequisites

• Python 3.x installed on local machine.
• A terminal, command prompt, or Jupyter Notebook environment.

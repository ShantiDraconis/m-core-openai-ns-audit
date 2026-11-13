# m-core-openai-ns-audit

O conteúdo verificável está em [`atk-kernel/`](atk-kernel/README.md): uma pasta Lean 4 independente (Mathlib `v4.34.1`) com lemas auxiliares do programa ATK sobre Navier–Stokes e auditoria automática de axiomas.

## Status

- `atk-kernel/` compila no CI (workflow `atk-kernel.yml`) e no Livro de Provas do MVEI (2026-11-11): 4 teoremas fechados, 0 condicionais, 0 abertos.
- `atk-kernel/scripts/check_axioms.py` só marca um teorema como CERTIFIED se o arquivo compila, os axiomas que ele usa estão contidos em {`propext`, `Classical.choice`, `Quot.sound`} (pode não usar nenhum), e a fonte não tem `sorry`, `admit` nem `axiom`. A lista esperada fica em `atk-kernel/scripts/expected_theorems.txt`.
- São lemas auxiliares; nenhum deles é a regularidade global de Navier–Stokes, que continua aberta.

`BACKWARD_ORACLE_PROTOCOL.md` descreve o protocolo de auditoria adotado nos repositórios do programa.

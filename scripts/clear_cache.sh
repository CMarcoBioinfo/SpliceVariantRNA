#!/bin/bash

read -t 21600 -p  "Voulez-vous vider la mémoire cache ? (y/n): " choice
if [ $? -ne 0 ]; then
    echo "Temps écoulé. Opération annulée : la mémoire cache n'a pas été vidée."
    exit 0
fi

if [[ "$choice" == "y" || "$choice" == "yes" || "$choice" == "o" || "$choice" == "oui" ]]; then
    if sudo -n true 2>/dev/null; then
        # Si l'utilisateur a déjà les droits sudo, vider la mémoire cache sans redemander le mot de passe
        echo "Vider la mémoire cache..."
        sudo -S sync && sudo -S sh -c 'echo 3 > /proc/sys/vm/drop_caches'
        echo "Mémoire cache vidée"
    else
        # Sinon, demander le mot de passe et vider la mémoire cache
        attempts=3
        while [ $attempts -gt 0 ]; do
            read -s -p "Entrer le mot de passe sudo pour vider la mémoire cache : " password
            echo
            echo $password | sudo -S sync && echo $password | sudo -S sh -c 'echo 3 > /proc/sys/vm/drop_caches'
            if [ $? -eq 0 ]; then
                echo "Mémoire cache vidée"
                exit 0
            else
                echo "Mot de passe incorrect. Il vous reste $((attempts-1)) tentatives."
            fi
            ((attempts--))
            sleep 1  # Ajout d'une pause pour donner du temps
        done
        echo "Trop de tentatives incorrectes. Opération annulée : la mémoire cache n'a pas été vidée."
    fi
else
    echo "Opération annulée : la mémoire cache n'a pas été vidée."
fi